# Covariance matrix construction for randomization inference
# Split from core_inference_functions.R (lines 80-269)
# Dependencies: stats (qwilcox), src/inference/statistics.R

#' Design-Based Covariance Bound Between Two Matched Sets (Theorem 1)
#'
#' Computes an upper bound on Cov(U_s, U_k) using the inverse-CDF / quantile
#' dot-product method. Uses a 10,000-point grid for numerical accuracy.
#'
#' @param nt_s Treated count in set s.
#' @param nc_s Control count in set s.
#' @param nt_k Treated count in set k.
#' @param nc_k Control count in set k.
#' @return Scalar covariance bound.
get_design_cov_bound <- function(nt_s, nc_s, nt_k, nc_k) {
  grid <- seq(0.0001, 0.9999, length.out = 10000)
  q_s  <- qwilcox(grid, m = nt_s, n = nc_s)
  q_k  <- qwilcox(grid, m = nt_k, n = nc_k)
  mean(q_s * q_k) - mean(q_s) * mean(q_k)
}


#' Build the Full Covariance Matrix for the Weighted U-Statistic
#'
#' Three covariance types:
#'   "naive"  — diagonal only (assumes independent matched sets)
#'   "design" — design-based worst-case off-diagonals (Theorem 1)
#'   "decay"  — design bound scaled by a network decay weighting matrix
#'
#' @param design_df  Data.frame with subclass and treatment columns.
#' @param weighting_mat S x S matrix of decay weights (required for "decay").
#' @param params  List with $type in c("naive", "design", "decay").
#' @param subclass_col Name of subclass column.
#' @param treat_col    Name of treatment column.
#' @return S x S covariance matrix Sigma.
get_U_covariance_matrix <- function(design_df, weighting_mat = NULL,
                                    params = list(type = "naive"),
                                    subclass_col = "subclass",
                                    treat_col = "Z") {
  valid_types <- c("naive", "design", "decay")
  if (is.null(params$type) || !(params$type %in% valid_types)) {
    stop(paste("params$type must be one of:", paste(valid_types, collapse = ", ")))
  }

  sc_vec <- design_df[[subclass_col]]
  z_vec  <- design_df[[treat_col]]
  sets   <- sort(unique(sc_vec))
  S      <- length(sets)

  set_stats <- lapply(sets, function(s) {
    idx <- which(sc_vec == s)
    nt  <- sum(z_vec[idx])
    nc  <- length(idx) - nt
    c(nt = nt, nc = nc)
  })

  Sigma <- matrix(0, S, S)

  for (i in seq_len(S)) {
    Sigma[i, i] <- get_U_moments(set_stats[[i]]["nt"], set_stats[[i]]["nc"])$var

    if (params$type != "naive" && i < S) {
      for (j in (i + 1):S) {
        cov_val <- get_design_cov_bound(
          set_stats[[i]]["nt"], set_stats[[i]]["nc"],
          set_stats[[j]]["nt"], set_stats[[j]]["nc"]
        )
        decay_factor <- if (params$type == "decay") weighting_mat[i, j] else 1
        Sigma[i, j] <- cov_val * decay_factor
        Sigma[j, i] <- cov_val * decay_factor
      }
    }
  }

  Sigma
}


#' Precompute the design-only covariance matrix (no decay weighting).
#'
#' Computes get_design_cov_bound for every pair of sets and diagonal variances.
#' The result can be reused across all (eta, rho) grid points via
#' get_U_cov_from_precomputed(), avoiding redundant computation.
#'
#' @param design_df  Data.frame with subclass and treatment columns.
#' @param subclass_col Name of subclass column.
#' @param treat_col    Name of treatment column.
#' @return S x S matrix where off-diagonals = cov_bound, diagonals = var.
precompute_design_cov_matrix <- function(design_df,
                                         subclass_col = "subclass",
                                         treat_col    = "Z") {
  sc_vec <- design_df[[subclass_col]]
  z_vec  <- design_df[[treat_col]]
  sets   <- sort(unique(sc_vec))
  S      <- length(sets)

  set_stats <- lapply(sets, function(s) {
    idx <- which(sc_vec == s)
    nt  <- sum(z_vec[idx])
    nc  <- length(idx) - nt
    c(nt = nt, nc = nc)
  })

  M <- matrix(0, S, S)
  for (i in seq_len(S)) {
    M[i, i] <- get_U_moments(set_stats[[i]]["nt"], set_stats[[i]]["nc"])$var
    if (i < S) {
      for (j in (i + 1):S) {
        v <- get_design_cov_bound(
          set_stats[[i]]["nt"], set_stats[[i]]["nc"],
          set_stats[[j]]["nt"], set_stats[[j]]["nc"]
        )
        M[i, j] <- v
        M[j, i] <- v
      }
    }
  }
  M
}


#' Build Sigma from a precomputed design matrix and a decay weighting matrix.
#'
#' Fast version: skips get_design_cov_bound entirely.
#' Sigma[i,j] = design_cov[i,j] * weighting_mat[i,j]  (off-diagonal)
#' Sigma[i,i] = design_cov[i,i]                         (diagonal unchanged)
#'
#' @param design_cov   S x S matrix from precompute_design_cov_matrix().
#' @param weighting_mat S x S decay weighting matrix.
#' @return S x S covariance matrix Sigma.
get_U_cov_from_precomputed <- function(design_cov, weighting_mat) {
  S     <- nrow(design_cov)
  diags <- diag(design_cov)
  Sigma <- design_cov * weighting_mat
  diag(Sigma) <- diags
  Sigma
}


#' Minimum-Distance Matrix Between Matched Sets (dmin method)
#'
#' D[s, k] = min distance between any unit in set s and any unit in set k,
#' using the unit-level distance matrix.
#'
#' @param unit_adj_mat n x n unit-level distance matrix (e.g. net_dist).
#' @param subclass_ids n-vector of subclass assignments.
#' @param method Currently only "dmin" is supported.
#' @return S x S distance matrix D with 0 on diagonal.
get_set_distance_matrix <- function(unit_adj_mat, subclass_ids, method = "dmin") {
  if (!method %in% c("dmin")) {
    stop("method must be 'dmin'")
  }

  sets     <- sort(unique(subclass_ids))
  S        <- length(sets)
  D        <- matrix(0, S, S)
  idx_list <- lapply(sets, function(s) which(subclass_ids == s))

  if (S < 2) return(D)

  for (i in seq_len(S - 1)) {
    for (j in (i + 1):S) {
      sub_vals <- unit_adj_mat[idx_list[[i]], idx_list[[j]], drop = FALSE]
      val <- suppressWarnings(min(sub_vals, na.rm = TRUE))
      if (!is.finite(val)) val <- Inf
      D[i, j] <- val
      D[j, i] <- val
    }
  }
  D
}


#' Exponential Decay Weighting Matrix Between Matched Sets (dmin-based)
#'
#' W[s, k] = eta * rho^(D[s,k] - 1)  if D[s,k] <= d0
#' W[s, k] = 0                         if D[s,k] > d0 or D[s,k] is Inf
#' W[s, s] = diag_val (default 1)
#'
#' @param dist_mat S x S set-distance matrix (from get_set_distance_matrix()).
#' @param eta  Decay strength at distance 1.
#' @param rho  Decay rate per unit distance.
#' @param d0   Max distance with nonzero weight (default Inf = no cutoff).
#' @param diag_val Value for diagonal entries.
#' @return S x S weighting matrix W.
get_set_decay_matrix <- function(dist_mat, eta, rho, d0 = Inf, diag_val = 1) {
  S <- nrow(dist_mat)
  W <- matrix(0, S, S)
  diag(W) <- diag_val

  if (S < 2) return(W)

  for (i in seq_len(S - 1)) {
    for (j in (i + 1):S) {
      d   <- dist_mat[i, j]
      if (!is.finite(d)) d <- Inf
      val <- eta * rho^(d - 1)
      if (!is.infinite(d0) && d > d0) val <- 0
      W[i, j] <- val
      W[j, i] <- val
    }
  }
  W
}


#' Mean Pairwise Decay Weighting Matrix (hard cutoff)
#'
#' W[s, k] = mean over all (i in s, j in k) of eta * rho^(d_ij - 1)
#' Pairs with d_ij > d0 contribute 0.
#'
#' @param unit_adj_mat n x n unit-level distance matrix.
#' @param subclass_ids n-vector of subclass assignments.
#' @param eta, rho, d0, diag_val As in get_set_decay_matrix().
#' @return S x S weighting matrix W.
get_set_mean_decay_matrix <- function(unit_adj_mat, subclass_ids,
                                       eta, rho, d0 = Inf, diag_val = 1) {
  sets     <- sort(unique(subclass_ids))
  S        <- length(sets)
  W        <- matrix(0, S, S)
  diag(W)  <- diag_val
  idx_list <- lapply(sets, function(s) which(subclass_ids == s))

  if (S < 2) return(W)

  for (i in seq_len(S - 1)) {
    for (j in (i + 1):S) {
      dists  <- unit_adj_mat[idx_list[[i]], idx_list[[j]], drop = FALSE]
      decays <- eta * rho^(dists - 1)
      if (!is.infinite(d0)) decays[dists > d0] <- 0
      val <- mean(decays, na.rm = TRUE)
      if (!is.finite(val)) val <- 0
      W[i, j] <- val
      W[j, i] <- val
    }
  }
  W
}


#' Mean Pairwise Decay Weighting Matrix (soft cutoff)
#'
#' Like get_set_mean_decay_matrix() but clamps distances at d0
#' instead of setting decays to 0: d_ij' = min(d_ij, d0).
#'
#' @inheritParams get_set_mean_decay_matrix
get_set_mean_decay_matrix_soft <- function(unit_adj_mat, subclass_ids,
                                            eta, rho, d0 = Inf, diag_val = 1) {
  sets     <- sort(unique(subclass_ids))
  S        <- length(sets)
  W        <- matrix(0, S, S)
  diag(W)  <- diag_val
  idx_list <- lapply(sets, function(s) which(subclass_ids == s))

  if (S < 2) return(W)

  for (i in seq_len(S - 1)) {
    for (j in (i + 1):S) {
      dists  <- unit_adj_mat[idx_list[[i]], idx_list[[j]], drop = FALSE]
      if (!is.infinite(d0)) dists <- pmin(dists, d0)
      decays <- eta * rho^(dists - 1)
      val    <- mean(decays, na.rm = TRUE)
      if (!is.finite(val)) val <- 0
      W[i, j] <- val
      W[j, i] <- val
    }
  }
  W
}
