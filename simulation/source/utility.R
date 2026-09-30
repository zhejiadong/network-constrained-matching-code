# Utility functions for matching
# Dependencies: dplyr, MASS, stats

#' Compute Mahalanobis Distances for Matching
#'
#' This function computes Mahalanobis distances between observations based on a
#' specified set of covariates. It can return either all-to-all distances or only
#' treated-to-control distances. The covariance matrix can be estimated either
#' using the overall sample covariance or a pooled within-group covariance
#' (treated and control groups separately, then combined).
#'
#' @param formula A formula of the form \code{Z ~ covariates}, where the left-hand
#'   side \code{Z} is a binary treatment indicator (0/1) and the right-hand side
#'   contains the covariates used to compute Mahalanobis distance. Factors are
#'   automatically dummy-encoded.
#' @param df A \code{data.frame} containing the treatment indicator and covariates.
#' @param mode Character string indicating which distances to return.
#'   \itemize{
#'     \item \code{"tc"} (default): return a treated × control distance matrix.
#'     \item \code{"all"}: return the full all-to-all distance matrix.
#'   }
#' @param cov.type Character string indicating how to estimate the covariance matrix.
#'   \itemize{
#'     \item \code{"pooled"} (default): pooled within-group covariance
#'       (treated and control separately, then combined).
#'     \item \code{"overall"}: covariance of the full sample (classical Mahalanobis).
#'   }
#' @param squared Logical; if \code{TRUE}, return squared Mahalanobis distances.
#'   If \code{FALSE} (default), return square-rooted distances (like
#'   \code{optmatch::match_on(method="mahalanobis")}).
#' @param id.col Optional string giving the column name in \code{df} to use as IDs
#'   for the rows/columns of the distance matrix. If \code{NULL} (default), row
#'   names from the model frame or observation indices are used.
#' @param standardize Default is true, in case some covariate has vary large variance,
#'   that dominate the distance measure
#' @return A numeric matrix of Mahalanobis distances. If \code{mode = "all"}, the
#'   matrix is square with all observations. If \code{mode = "tc"}, the rows are
#'   treated units and the columns are control units.
maha_my <- 
  function(formula, df, mode = c("tc","all"), cov.type = c("pooled","overall"),
           squared = FALSE, id.col = NULL, standardize = TRUE){
  
  mode <- match.arg(mode) 
  cov.type <- match.arg(cov.type)
  
  mf <- model.frame(formula, data = df, na.action = na.omit)
  Z <- model.response(mf)
  rhs_terms <- reformulate(attr(terms(formula), "term.labels"))
  #model matrix prepare
  X <- model.matrix(rhs_terms, data = mf)
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
  }
  if (standardize) {
    X <- scale(X, center = TRUE, scale = TRUE)
  } else {
    X <- scale(X, center = TRUE, scale = FALSE)
  }
  
  ### Ids
  ids <- rownames(mf)
  if (is.null(ids)) ids <- seq_len(nrow(mf))
  rownames(X) <- ids
  
  ### calculate covariance, if pooled or overall.
  if (cov.type == "pooled") {
    n1 <- sum(Z == 1); n0 <- sum(Z == 0)
    S1 <- stats::cov(X[Z == 1, , drop = FALSE])
    S0 <- stats::cov(X[Z == 0, , drop = FALSE])
    S <- ((n1 - 1) * S1 + (n0 - 1) * S0) / (n1 + n0 - 2)
  } else if (cov.type == "overall") {
    S <- stats::cov(X)
  }
  
  Sinv <- MASS::ginv(S)
  Q  <- X %*% Sinv %*% t(X)
  d2 <- outer(diag(Q), diag(Q), "+") - 2 * Q
  d2[d2 < 0] <- 0
  
  D <- if (squared) d2 else sqrt(d2)
  
  if (mode == "all") {
    rownames(D) <- colnames(D) <- ids
    return(D)
  } else {
    treat_ids <- ids[Z == 1]
    ctrl_ids  <- ids[Z == 0]
    Dtc <- D[ids %in% treat_ids, ids %in% ctrl_ids, drop = FALSE]
    rownames(Dtc) <- treat_ids
    colnames(Dtc) <- ctrl_ids
    return(Dtc)
  }
}


#' Apply network-based penalty to a distance matrix
#'
#' This function adjusts a Mahalanobis distance matrix by applying
#' a penalty based on network proximity. Units that are too distant
#' in the network (according to a specified prox_thresh) will have their
#' distance set to a large penalty value.
#'
#' @param maha_dist A numeric distance matrix (e.g., from \code{maha_my()}).
#' @param net_dist A network distance matrix, typically obtained from
#'   \code{igraph::distances()}.
#' @param prox_thresh A numeric threshold for network proximity. Pairs with
#'   proximity AT OR ABOVE this value (i.e., close in the network) are penalized.
#'   Default is 0.25.
#' @param penalty A numeric penalty value assigned to disallowed matches.
#'   If \code{NULL}, the penalty is set to \code{max(maha_dist, na.rm=TRUE) * 1000}.
#'    Do not use Inf in case of unexpected error when using quickmatch, and so
#'    may in matchit/optmatch
#' @param treat treatment assignment variable
#' @return A numeric matrix of the same size as \code{maha_dist}, with penalties applied.
#'   Row and column names are preserved.
network_penalty <- function(maha_dist, net_dist, prox_thresh = 0.25, lessequal = FALSE,
                            penalty = NULL, treat) {
  
  # 1. Compute network proximity (inverse of network distance)
  net_prox <- 1 / net_dist
  # 2. Define penalty if not provided
  if (is.null(penalty)) {
    penalty <- max(maha_dist[is.finite(maha_dist)], na.rm = TRUE) * 100000000
  }
  
  if(nrow(maha_dist) != ncol(maha_dist)){
    net_prox <- net_prox[ which(treat == 1), which(treat == 0)]
  }else{
    net_prox =   net_prox
  }

    # 3. Apply penalty rule
    if (lessequal) {
      maha_penalty <- ifelse(net_prox <= prox_thresh, maha_dist, penalty)
    }else{
      maha_penalty <- ifelse(net_prox < prox_thresh, maha_dist, penalty)
    }

  
 
  
  # 4. Preserve dimnames
  rownames(maha_penalty) <- rownames(maha_dist)
  colnames(maha_penalty) <- colnames(maha_dist)
  
  return(maha_penalty)
}