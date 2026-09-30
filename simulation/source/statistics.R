# U-statistic computation for randomization inference
# Split from core_inference_functions.R (lines 1-78)
# Dependencies: base R only

#' Calculate Observed Weighted U-Statistics Across Matched Sets
#'
#' Computes the Mann-Whitney U-statistic within each matched set, then
#' aggregates to a single weighted test statistic T_obs.
#'
#' @param data Data.frame with outcome, treatment, and subclass columns.
#' @param weight_type Weighting scheme:
#'   "ns"  → w_s = 1/(n_s + 1)  (default, Wilcoxon-type)
#'   "ntc" → w_s = 1/(n_t * n_c)
#' @param subclass_col Name of the subclass/stratum column.
#' @param treat_col    Name of the binary treatment column (0/1).
#' @param outcome_col  Name of the outcome column.
#'
#' @return A list with:
#'   - set_stats: data.frame (one row per set) with subclass, n, nt, nc,
#'                U_obs, E_null, weight
#'   - T_obs: scalar weighted sum of observed U-statistics
#'   - E_T:   scalar weighted sum of null expectations
get_observed_stats <- function(data, weight_type = "ns",
                               subclass_col = "subclass",
                               treat_col = "Z",
                               outcome_col = "Y") {
  data <- data[order(data[[subclass_col]]), ]
  sets <- sort(unique(data[[subclass_col]]))

  results_list <- lapply(sets, function(s) {
    sub <- data[data[[subclass_col]] == s, ]
    Y  <- sub[[outcome_col]]
    Z  <- sub[[treat_col]]
    n  <- nrow(sub)
    nt <- sum(Z)
    nc <- n - nt

    U_val <- if (nt == 0 || nc == 0) {
      0
    } else {
      sum(outer(Y[Z == 1], Y[Z == 0], function(t, c) as.integer(t > c)))
    }

    E_null <- (nt * nc) / 2

    w <- switch(weight_type,
      ns  = 1 / (n + 1),
      ntc = if (nt > 0 && nc > 0) 1 / (nt * nc) else 0,
      stop("weight_type must be 'ns' or 'ntc'")
    )

    data.frame(subclass = s, n = n, nt = nt, nc = nc,
               U_obs = U_val, E_null = E_null, weight = w)
  })

  results_df <- do.call(rbind, results_list)
  T_obs <- sum(results_df$weight * results_df$U_obs)
  E_T   <- sum(results_df$weight * results_df$E_null)

  list(set_stats = results_df, T_obs = T_obs, E_T = E_T)
}


#' Null Mean and Variance of the Mann-Whitney U-Statistic
#'
#' Under the sharp null of no treatment effect, within one matched set:
#'   E[U] = nt * nc / 2
#'   Var[U] = nt * nc * (nt + nc + 1) / 12
#'
#' @param nt Number of treated units.
#' @param nc Number of control units.
#' @return List with $mu (mean) and $var (variance).
get_U_moments <- function(nt, nc) {
  list(
    mu  = (nt * nc) / 2,
    var = (nt * nc * (nt + nc + 1)) / 12
  )
}
