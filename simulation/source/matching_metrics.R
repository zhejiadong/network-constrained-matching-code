# Matching evaluation metric functions
#
# All functions are pure (no I/O, no global state).
# Dependencies: MASS (for ginv fallback)

# =============================================================================
# 1. Covariate Balance (SMD)
# =============================================================================

#' Unmatched (before) SMD for one covariate
#'
#' Uses pooled SD: sqrt((var_t + var_c) / 2)
#'
#' @param x Numeric covariate vector.
#' @param z Binary treatment vector (0/1).
#' @return Scalar signed SMD (use abs() for |SMD|).
smd_before <- function(x, z) {
  xt <- x[z == 1]; xc <- x[z == 0]
  if (length(xt) < 2 || length(xc) < 2) return(NA_real_)
  s <- sqrt((var(xt) + var(xc)) / 2)
  if (!is.finite(s) || s == 0) return(NA_real_)
  (mean(xt) - mean(xc)) / s
}


#' Matched (after) SMD for one covariate
#'
#' Computes within-set mean differences, aggregates by set size,
#' standardizes by the global pooled SD.
#'
#' @param df Data.frame with covariate, subclass, and treatment columns.
#' @param xcol Covariate column name.
#' @param scol Subclass column name.
#' @param zcol Treatment column name.
#' @return Scalar signed matched SMD.
smd_after <- function(df, xcol, scol, zcol = "Z") {
  x <- df[[xcol]]; z <- df[[zcol]]; s <- df[[scol]]
  sd0 <- sqrt((var(x[z == 1], na.rm = TRUE) + var(x[z == 0], na.rm = TRUE)) / 2)
  if (!is.finite(sd0) || sd0 == 0) return(NA_real_)

  sets <- split(seq_len(nrow(df)), s)
  num <- 0; den <- 0
  for (ids in sets) {
    zz <- z[ids]
    if (!any(zz == 1) || !any(zz == 0)) next
    d   <- mean(x[ids][zz == 1], na.rm = TRUE) - mean(x[ids][zz == 0], na.rm = TRUE)
    w   <- length(ids)
    num <- num + w * d
    den <- den + w
  }
  if (den == 0) return(NA_real_)
  (num / den) / sd0
}


#' Compute |SMD| for all covariates, before and after matching
#'
#' @param df Data.frame with covariates, treatment, and subclass columns.
#' @param covars Character vector of covariate names.
#' @param subclass_methods Named list: method_label → subclass_column.
#' @param zcol Treatment column name.
#' @return Long data.frame with columns: method, covariate, abs_smd.
compute_smd_long <- function(df, covars, subclass_methods, zcol = "Z") {
  rows <- list()
  k <- 1L

  for (x in covars) {
    rows[[k]] <- data.frame(method = "Before", covariate = x,
                             abs_smd = abs(smd_before(df[[x]], df[[zcol]])),
                             stringsAsFactors = FALSE)
    k <- k + 1L
  }

  for (mname in names(subclass_methods)) {
    scol <- subclass_methods[[mname]]
    if (!scol %in% names(df) || all(is.na(df[[scol]]))) next
    for (x in covars) {
      rows[[k]] <- data.frame(method = mname, covariate = x,
                               abs_smd = abs(smd_after(df, x, scol, zcol)),
                               stringsAsFactors = FALSE)
      k <- k + 1L
    }
  }

  do.call(rbind, rows)
}


# =============================================================================
# 2. Within-Set Mahalanobis Distance (covariate similarity)
# =============================================================================

#' Compute global inverse covariance matrix (with ginv fallback)
.get_S_inv <- function(X) {
  S <- cov(X)
  if (abs(det(S)) > 1e-10) solve(S) else MASS::ginv(S)
}


#' Average Mahalanobis distance between all treated-control pairs (before matching)
#'
#' @param df Data.frame with covariate and treatment columns.
#' @param covars Character vector of covariate names.
#' @param zcol Treatment column name.
#' @return Scalar mean Mahalanobis distance.
mahal_before <- function(df, covars, zcol = "Z") {
  X     <- as.matrix(df[, covars, drop = FALSE])
  S_inv <- .get_S_inv(X)
  z     <- df[[zcol]]
  tr    <- which(z == 1)
  co    <- which(z == 0)

  dists <- vapply(tr, function(ti) {
    mean(vapply(co, function(ci) {
      diff <- X[ti, ] - X[ci, ]
      sqrt(as.numeric(t(diff) %*% S_inv %*% diff))
    }, numeric(1)))
  }, numeric(1))

  mean(dists, na.rm = TRUE)
}


#' Average within-set Mahalanobis distance between treated-control pairs
#'
#' Uses global pooled covariance for standardization.
#'
#' @param df Data.frame with covariate, treatment, and subclass columns.
#' @param scol Subclass column name.
#' @param covars Character vector of covariate names.
#' @param zcol Treatment column name.
#' @return Scalar mean within-set Mahalanobis distance.
mahal_after <- function(df, scol, covars, zcol = "Z") {
  X     <- as.matrix(df[, covars, drop = FALSE])
  S_inv <- .get_S_inv(X)
  z     <- df[[zcol]]
  sets  <- split(seq_len(nrow(df)), df[[scol]])

  pair_dists <- unlist(lapply(sets, function(ids) {
    zz <- z[ids]
    tr <- ids[zz == 1]
    co <- ids[zz == 0]
    if (length(tr) == 0 || length(co) == 0) return(numeric(0))
    vapply(tr, function(ti) {
      vapply(co, function(ci) {
        diff <- X[ti, ] - X[ci, ]
        sqrt(as.numeric(t(diff) %*% S_inv %*% diff))
      }, numeric(1))
    }, numeric(length(co)))
  }))

  mean(pair_dists, na.rm = TRUE)
}


#' Compute within-set Mahalanobis for all methods
#'
#' @param df Data.frame.
#' @param covars Covariate names.
#' @param subclass_methods Named list: label → subclass column.
#' @param zcol Treatment column.
#' @return Long data.frame with columns: method, mahal_dist.
compute_mahal_long <- function(df, covars, subclass_methods, zcol = "Z") {
  rows <- list(
    data.frame(method = "Before", mahal_dist = mahal_before(df, covars, zcol),
               stringsAsFactors = FALSE)
  )
  for (mname in names(subclass_methods)) {
    scol <- subclass_methods[[mname]]
    if (!scol %in% names(df) || all(is.na(df[[scol]]))) next
    rows[[length(rows) + 1L]] <- data.frame(
      method = mname,
      mahal_dist = mahal_after(df, scol, covars, zcol),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}


# =============================================================================
# 3. Network Distance Statistics
# =============================================================================

#' Pairwise distance stats within or across matched sets
#'
#' @param df Data.frame with subclass column.
#' @param scol Subclass column name.
#' @param dmat n x n network distance matrix.
#' @param scope "within" (same set) or "across" (different sets).
#' @return Named numeric: avg, p1 (d=1), p2 (d=2), pge3 (d>=3).
net_pair_stats <- function(df, scol, dmat, scope = c("within", "across"),
                           max_across = 5000L) {
  scope <- match.arg(scope)
  s <- df[[scol]]
  pr <- combn(seq_len(nrow(df)), 2)
  same <- s[pr[1, ]] == s[pr[2, ]]
  keep <- if (scope == "within") same else !same
  pr <- pr[, keep, drop = FALSE]

  if (ncol(pr) == 0) return(c(avg = NA_real_, p1 = NA_real_, p2 = NA_real_, pge3 = NA_real_))

  # For across pairs, subsample to avoid excessive computation
  if (scope == "across" && ncol(pr) > max_across) {
    idx <- sample.int(ncol(pr), max_across)
    pr  <- pr[, idx, drop = FALSE]
  }

  d_all <- dmat[cbind(pr[1, ], pr[2, ])]
  d_all <- d_all[is.finite(d_all)]
  if (!length(d_all)) return(c(avg = NA_real_, p1 = NA_real_, p2 = NA_real_, pge3 = NA_real_))

  c(avg = mean(d_all), p1 = mean(d_all == 1), p2 = mean(d_all == 2), pge3 = mean(d_all >= 3))
}


#' Minimum-distance stats across matched sets
#'
#' For each pair of matched sets (A, B), compute
#'   d_min(A, B) = min_{i in A, j in B} d(i, j)
#' and summarize the distribution of these set-pair minima.
#'
#' @param df Data.frame with subclass column.
#' @param scol Subclass column name.
#' @param dmat n x n network distance matrix.
#' @return Named numeric: avg, p1 (d=1), p2 (d=2), pge3 (d>=3).
net_across_dmin_stats <- function(df, scol, dmat) {
  df_sub <- df[!is.na(df[[scol]]), , drop = FALSE]
  if (nrow(df_sub) == 0) {
    return(c(avg = NA_real_, p1 = NA_real_, p2 = NA_real_, pge3 = NA_real_))
  }

  df_sub$subclass_tmp <- df_sub[[scol]]
  sets <- sort(unique(df_sub$subclass_tmp))
  if (length(sets) < 2) {
    return(c(avg = NA_real_, p1 = NA_real_, p2 = NA_real_, pge3 = NA_real_))
  }

  get_idx <- function(sub_df) as.integer(rownames(sub_df))
  dmins <- c()

  for (i in seq_along(sets)[-length(sets)]) {
    for (j in (i + 1):length(sets)) {
      sub_a <- df_sub[df_sub$subclass_tmp == sets[i], , drop = FALSE]
      sub_b <- df_sub[df_sub$subclass_tmp == sets[j], , drop = FALSE]
      idx_a <- get_idx(sub_a)
      idx_b <- get_idx(sub_b)

      valid <- all(idx_a %in% seq_len(nrow(dmat))) &&
               all(idx_b %in% seq_len(nrow(dmat)))
      if (!valid || length(idx_a) == 0 || length(idx_b) == 0) next

      dmins <- c(dmins, min(dmat[idx_a, idx_b, drop = FALSE], na.rm = TRUE))
    }
  }

  dmins <- dmins[is.finite(dmins)]
  if (!length(dmins)) {
    return(c(avg = NA_real_, p1 = NA_real_, p2 = NA_real_, pge3 = NA_real_))
  }

  c(avg = mean(dmins), p1 = mean(dmins == 1), p2 = mean(dmins == 2), pge3 = mean(dmins >= 3))
}


#' All-pairs distance stats (before matching)
#'
#' @param df Data.frame.
#' @param dmat n x n network distance matrix.
#' @return Named numeric: avg, p1, p2, pge3.
net_before_stats <- function(df, dmat) {
  pr <- combn(seq_len(nrow(df)), 2)
  d  <- dmat[cbind(pr[1, ], pr[2, ])]
  d  <- d[is.finite(d)]
  c(avg = mean(d), p1 = mean(d == 1), p2 = mean(d == 2), pge3 = mean(d >= 3))
}


#' Compute network distance stats for all methods
#'
#' @param df Data.frame.
#' @param dmat n x n network distance matrix.
#' @param subclass_methods Named list: label → subclass column.
#' @return Long data.frame with columns: method, scope, avg, p1, p2, pge3.
compute_net_long <- function(df, dmat, subclass_methods) {
  b <- net_before_stats(df, dmat)
  rows <- list(data.frame(method = "Before", scope = "all",
                           avg = b["avg"], p1 = b["p1"], p2 = b["p2"], pge3 = b["pge3"],
                           stringsAsFactors = FALSE))

  for (mname in names(subclass_methods)) {
    scol <- subclass_methods[[mname]]
    if (!scol %in% names(df) || all(is.na(df[[scol]]))) next

    wi <- net_pair_stats(df, scol, dmat, "within")
    ac <- net_pair_stats(df, scol, dmat, "across")
    dm <- net_across_dmin_stats(df, scol, dmat)

    rows[[length(rows) + 1L]] <- data.frame(
      method = mname, scope = "within",
      avg = wi["avg"], p1 = wi["p1"], p2 = wi["p2"], pge3 = wi["pge3"],
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1L]] <- data.frame(
      method = mname, scope = "across",
      avg = ac["avg"], p1 = ac["p1"], p2 = ac["p2"], pge3 = ac["pge3"],
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1L]] <- data.frame(
      method = mname, scope = "across_dmin",
      avg = dm["avg"], p1 = dm["p1"], p2 = dm["p2"], pge3 = dm["pge3"],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}
