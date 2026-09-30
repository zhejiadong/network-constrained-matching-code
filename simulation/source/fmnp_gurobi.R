# Full Matching with Network Penalty using Gurobi optimizer
# Dependencies: Matrix, gurobi, igraph
library(Matrix)
library(gurobi)
make_pairs <- function(keep_slice, prox_mat, prox_thresh) {
  cand <- which(keep_slice)  
  if (length(cand) < 2) return(data.frame(c1 = integer(0), c2 = integer(0)))
  prox_sub <- prox_mat[cand, cand, drop = FALSE]
  diag(prox_sub) <- 0
  close_ij <- which(prox_sub >= prox_thresh & upper.tri(prox_sub), arr.ind = TRUE)
  
  if (nrow(close_ij) == 0) {
    return(data.frame(c1 = integer(0), c2 = integer(0)))
  }
  data.frame(
    c1 = cand[close_ij[, "row"]],
    c2 = cand[close_ij[, "col"]]
  )
}

extract_subclasses_edges <- function(res, t_ids, c_ids, E){
  nT <- length(t_ids); nC <- length(c_ids)
  sol <- res$x
  chosen <- which(sol > 0.5)
  sel_edges <- E[chosen, ]
  edges <- as.vector(t(cbind(sel_edges$t, nT + sel_edges$c)))
  g.sub <- igraph::make_empty_graph(n = nT + nC, directed = FALSE)
  g.sub <- igraph::add_edges(g.sub, edges)
  memb <- igraph::components(g.sub)$membership
  remap <- match(memb, unique(memb))
  
  N <- max(c(t_ids, c_ids))
  subclass <- rep(NA_integer_, N)
  subclass[t_ids] <- remap[1:nT]
  subclass[c_ids] <- remap[(nT + 1):(nT + nC)]
  factor(subclass)
}

build_gurobi_model <- function(E, eid,
                               nT, nC,
                               minC_per_T, maxC_per_T,
                               minT_per_C, maxT_per_C,
                               cc_pairs_list, tt_pairs_list) {
  
  m <- nrow(E)            # number of edge variables
  obj <- E$cost
  vtype <- rep('B', m)    # binary vars
  
  # We'll accumulate a sparse constraint matrix in triplet form
  Ai <- integer(0); Aj <- integer(0); Ax <- numeric(0)
  sense <- character(0); rhs <- numeric(0)
  
  add_row <- function(cols, coeffs, s, r) {
    # cols are variable indices (1..m), coeffs matching numeric coefficients
    # s is one of c('<=','>=','='), r is numeric rhs
    row_idx <- length(rhs) + 1L
    Ai <<- c(Ai, rep(row_idx, length(cols)))
    Aj <<- c(Aj, cols)
    Ax <<- c(Ax, coeffs)
    sense <<- c(sense, s)
    rhs <<- c(rhs, r)
  }
  
  ## 1) Per-treated lower/upper bounds: sum_e x_e >= minC_per_T  and <= maxC_per_T
  for (t in 1:nT) {
    cols <- E$e[E$t == t]
    if (length(cols)) {
      if (is.finite(minC_per_T) && minC_per_T > 0) {
        add_row(cols, rep(1, length(cols)), ">=", minC_per_T)
      }
      if (is.finite(maxC_per_T)) {
        add_row(cols, rep(1, length(cols)), "<=", maxC_per_T)
      }
    } else {
      # If a treated has no feasible edges, the model is infeasible.
    }
  }
  
  ## 2) Per-control lower/upper bounds: sum_e x_e >= minT_per_C and <= maxT_per_C
  for (c in 1:nC) {
    cols <- E$e[E$c == c]
    if (length(cols)) {
      if (is.finite(minT_per_C) && minT_per_C > 0) {
        add_row(cols, rep(1, length(cols)), ">=", minT_per_C)
      }
      if (is.finite(maxT_per_C)) {
        add_row(cols, rep(1, length(cols)), "<=", maxT_per_C)
      }
    }
  }
  
  cc_pairs_mat_list <- lapply(seq_len(nT), function(t) {
    cp <- cc_pairs_list[[t]]
    if (is.null(cp) || nrow(cp) == 0) return(NULL)
    e1 <- eid[t, cp$c1]
    e2 <- eid[t, cp$c2]
    keep <- (e1 > 0L) & (e2 > 0L)
    if (!any(keep)) return(NULL)
    cbind(e1[keep], e2[keep])  # each row = (e1, e2)
  })
  cc_pairs_mat <- if (length(cc_pairs_mat_list)) do.call(rbind, cc_pairs_mat_list) else NULL
  
  if (!is.null(cc_pairs_mat)) {
    k <- nrow(cc_pairs_mat)
    new_rows <- (length(rhs) + 1L):(length(rhs) + k)
    Ai <- c(Ai, rep(new_rows, each = 2L))
    Aj <- c(Aj, as.vector(t(cc_pairs_mat)))
    Ax <- c(Ax, rep(1, 2L * k))
    sense <- c(sense, rep("<=", k))
    rhs <- c(rhs, rep(1, k))
  }
  
  ## Build all T–T constraints at once
  tt_pairs_mat_list <- lapply(seq_len(nC), function(c) {
    tp <- tt_pairs_list[[c]]
    if (is.null(tp) || nrow(tp) == 0) return(NULL)
    e1 <- eid[tp$c1, c]
    e2 <- eid[tp$c2, c]
    keep <- (e1 > 0L) & (e2 > 0L)
    if (!any(keep)) return(NULL)
    cbind(e1[keep], e2[keep])  # each row = (e1, e2)
  })
  tt_pairs_mat <- if (length(tt_pairs_mat_list)) do.call(rbind, tt_pairs_mat_list) else NULL
  
  if (!is.null(tt_pairs_mat)) {
    k <- nrow(tt_pairs_mat)
    new_rows <- (length(rhs) + 1L):(length(rhs) + k)
    Ai <- c(Ai, rep(new_rows, each = 2L))
    Aj <- c(Aj, as.vector(t(tt_pairs_mat)))
    Ax <- c(Ax, rep(1, 2L * k))
    sense <- c(sense, rep("<=", k))
    rhs <- c(rhs, rep(1, k))
  }
  
  # Build sparse A
  if (length(Ai) == 0L) {
    A <- Matrix::Matrix(0, nrow = 0, ncol = m, sparse = TRUE)
  } else {
    A <- Matrix::sparseMatrix(i = Ai, j = Aj, x = Ax, dims = c(length(rhs), m))
  }
  
  model <- list(
    A = A,
    obj = obj,
    sense = sense,
    rhs = rhs,
    vtype = vtype,
    modelsense = "min"
  )
  
  return(model)
}


#' Full Matching with Network Penalty (Gurobi)
#'
#' @param formula A formula of the form Z ~ covariates.
#' @param data Data frame containing treatment, covariates, and network IDs.
#' @param D Distance matrix (treated × control) after network penalty.
#' @param net_proximity Full network proximity matrix (inverse of distances).
#' @param caliper Caliper on Mahalanobis distance.
#' @param prox_thresh Proximity threshold for network penalty.
#' @param min.controls Minimum number of controls per treated (can be fractional).
#' @param max.controls Maximum number of controls per treated.
#' @param timelimit Time limit for Gurobi (seconds).
#' @param mipgap MIP gap tolerance.
#' @param thread Number of threads for Gurobi.
#' @param OutputFlag Gurobi output flag (0 = silent).
#'
#' @return A list with components:
#'   \item{data}{Original data frame with an added 'subclass' factor column.}
#'   \item{res}{Raw Gurobi solution object.}
FMNP_gurobi <- function(formula, data, D, net_proximity,
                        caliper = 4.5, prox_thresh = 0.5,
                        min.controls = 0.1, max.controls = 50,
                        timelimit = 600, mipgap = 0.01, thread = 1L,
                        OutputFlag = 0) {
  ## step 1: setup basic
  mf <- model.frame(formula, data)
  Z <- model.response(mf)
  if (!all(Z %in% c(0,1))) stop("Treatment must be 0/1.")
  t_ids <- which(Z == 1)
  c_ids <- which(Z == 0)
  nT <- length(t_ids)
  nC <- length(c_ids)
  if (nT == 0 || nC == 0) stop("Need at least one treated and one control.")
  
  prox_cc <- as.matrix(net_proximity[c_ids, c_ids, drop = FALSE])
  prox_tc <- as.matrix(net_proximity[t_ids, c_ids, drop = FALSE])
  prox_tt <- as.matrix(net_proximity[t_ids, t_ids, drop = FALSE])
  keep <- (D < caliper) & (prox_tc < prox_thresh)
  if (!any(keep)) stop("No feasible T–C edges under caliper/proximity.")
  
  deg_T <- rowSums(keep); deg_C <- colSums(keep)
  if (any(deg_T == 0)) stop("Infeasible: some treated have no feasible control (keep row = 0).")
  # If we allow controls to be unused (min.controls == 0), deg_C == 0 is OK.
  if (min.controls != 0 && any(deg_C == 0)) {
    stop("Infeasible: some controls have no feasible treated (keep col = 0).")
  }
  
  # number of control per treated
  # NOTE: In full matching / FMNP we still require each treated to have at least 1 control.
  minC_per_T <- max(1, min.controls)
  maxC_per_T <- max.controls
  
  # number of treated per control
  # If min.controls == 0, allow controls to be unused (minT_per_C = 0) and remove the upper bound.
  if (min.controls == 0) {
    minT_per_C <- 0
    maxT_per_C <- Inf
  } else {
    minT_per_C <- 1
    maxT_per_C <- floor(1 / min.controls)
    if (!is.finite(maxT_per_C)) maxT_per_C <- nT
  }
  cc_pairs_list <- lapply(1:nT, function(t) make_pairs(keep[t, ], prox_cc, prox_thresh))
  tt_pairs_list <- lapply(1:nC, function(c) make_pairs(keep[, c], prox_tt, prox_thresh))
  
  ## step 2: build edge list
  which_keep <- which(keep, arr.ind = TRUE)
  E <- data.frame(
    e = seq_len(nrow(which_keep)),
    t = which_keep[, 1],
    c = which_keep[, 2]
  )
  E$cost <- D[cbind(E$t, E$c)]
  eid <- matrix(0L, nrow = nT, ncol = nC)
  eid[cbind(E$t, E$c)] <- E$e
  nE <- nrow(E)
  
  ## step 3: gurobi model
  model <- build_gurobi_model(E, eid, nT, nC, 
                              minC_per_T, maxC_per_T,
                              minT_per_C, maxT_per_C,
                              cc_pairs_list, tt_pairs_list)
  
  res <- gurobi::gurobi(
    model,
    params = list(
      TimeLimit = timelimit,
      MIPGap = mipgap,
      Threads = as.integer(max(1, thread)),
      OutputFlag = OutputFlag
    )
  )

  ## step 4: check solution status
  if (res$status == "INFEASIBLE") stop("Gurobi: model is infeasible.")
  if (res$status == "INF_OR_UNBD") stop("Gurobi: model is infeasible or unbounded.")
  if (!res$status %in% c("OPTIMAL", "TIME_LIMIT") || is.null(res$x)) {
    stop(sprintf("Gurobi returned status '%s' with no solution.", res$status))
  }

  ## step 5: extract solution
  data$subclass <- extract_subclasses_edges(res, t_ids, c_ids, E)
  return(list(data = data, res = res))
}