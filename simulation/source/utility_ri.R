# Main function for randomization inference with normal approximation
# Dependencies: dplyr, stats, statistics.R, covariance.R

#' Run Randomization Inference
#' 
#' @param df Dataframe containing outcome, treatment, and subclass.
#' @param net_dist N x N unit-level distance matrix (Required if cov_type = "decay").
#' @param cov_type Type of dependence: "naive", "design", or "decay".
#' @param decay_params List containing 'eta', 'rho', 'd0', optional 'set_dist' in {"dmin","mean","mean_soft"}.
#' @param weight_type "ns" or "ntc" (default "ns").
#' @param subclass_col Name of subclass column.
#' @param treat_col Name of treatment column.
#' @param outcome_col Name of outcome column.
#' 
#' @return A list containing:
#'   - summary: Dataframe with T_obs, Expectation, Variance, Z-score, P-value.
#'   - Sigma: The calculated covariance matrix used.
#'   - stats_detail: Detailed U-statistics per set.
run_randomization_inference <- function(df, 
                                        net_dist = NULL,
                                        cov_type = c("naive", "design", "decay"),
                                        decay_params = list(eta = 1, rho = 1, d0 = Inf, set_dist = "dmin"),
                                        weight_type = "ns",
                                        subclass_col = "subclass",
                                        treat_col = "Z",
                                        outcome_col = "Y") {
  
  cov_type <- match.arg(cov_type)
  # --- Step 1: Calculate Observed Statistics (Script 2) ---
  # This gives us T_obs, E[T], and the weights vector
  obs_res <- get_observed_stats(
    data = df, 
    weight_type = weight_type,
    subclass_col = subclass_col, 
    treat_col = treat_col,        # Note: Script 2 uses 'trt_col'
    outcome_col = outcome_col
  )
  #dist_mat_sets <- NULL
  # Extract weights vector (ensure order matches covariance matrix)
  weights_vec <- obs_res$set_stats$weight
  
  
  
  # --- Step 2: Calculate Covariance Matrix (Script 1) ---
  cov_params <- list(type = cov_type)
  weighting_matrix <- NULL

  if (cov_type == "decay") {
    if (is.null(net_dist)) stop("Error: 'net_dist' (N x N) required for decay.")
    if (is.null(decay_params)) stop("Error: no decay parameters; please specify list(eta, rho, d0, set_dist)")

    if (is.null(decay_params$eta)) decay_params$eta <- 1
    if (is.null(decay_params$rho)) decay_params$rho <- 1
    if (is.null(decay_params$d0))  decay_params$d0  <- Inf
    if (is.null(decay_params$set_dist)) decay_params$set_dist <- "dmin"

    if (!(decay_params$set_dist %in% c("dmin", "mean", "mean_soft"))) {
      stop("decay_params$set_dist must be one of: 'dmin', 'mean', 'mean_soft'")
    }

    # Decay: Sigma_sk = design_bound_sk * weight_sk
    if (decay_params$set_dist == "mean") {
      weighting_matrix <- get_set_mean_decay_matrix(
        unit_adj_mat = net_dist,
        subclass_ids = df[[subclass_col]],
        eta = decay_params$eta,
        rho = decay_params$rho,
        d0  = decay_params$d0,
        diag_val = 1
      )
    } else if (decay_params$set_dist == "mean_soft") {
      weighting_matrix <- get_set_mean_decay_matrix_soft(
        unit_adj_mat = net_dist,
        subclass_ids = df[[subclass_col]],
        eta = decay_params$eta,
        rho = decay_params$rho,
        d0  = decay_params$d0,
        diag_val = 1
      )
    } else {
      set_dist_mat <- get_set_distance_matrix(
        unit_adj_mat = net_dist,
        subclass_ids = df[[subclass_col]],
        method = "dmin"
      )
      weighting_matrix <- get_set_decay_matrix(
        dist_mat = set_dist_mat,
        eta = decay_params$eta,
        rho = decay_params$rho,
        d0  = decay_params$d0,
        diag_val = 1
      )
    }

    Sigma <- get_U_covariance_matrix(
      design_df = df,
      weighting_mat = weighting_matrix,
      params = cov_params,
      subclass_col = subclass_col,
      treat_col = treat_col
    )
  } else {
    Sigma <- get_U_covariance_matrix(
      design_df = df,
      weighting_mat = weighting_matrix,
      params = cov_params,
      subclass_col = subclass_col,
      treat_col = treat_col
    )
  }
  # --- Step 3: Inference (Normal Approximation) ---
  
  # Variance of Weighted Sum: w' * Sigma * w
  Var_T <- as.numeric(t(weights_vec) %*% Sigma %*% weights_vec)
  
  # Guard against negative variance (can occur with aggressive decay parameters)
  if (Var_T < 0) {
    stop(sprintf(
      "Calculated variance is negative (Var_T = %.6g). ",
      "Check decay parameters (eta, rho, d0) or the distance matrix.",
      Var_T
    ))
  } else {
    # Z-Score
    Z_score <- (obs_res$T_obs - obs_res$E_T) / sqrt(Var_T)
    
    # Two-Sided P-Value
    p_val <- 2 * (1 - pnorm(abs(Z_score)))
  }
  
  # --- Step 4: Formatting Results ---
  summary_df <- data.frame(
    Method      = cov_type,
    Decay_Set_Dist = if (cov_type == "decay") decay_params$set_dist else "N/A",
    T_observed  = obs_res$T_obs,
    Expectation = obs_res$E_T,
    Variance    = Var_T,
    Std_Error   = sqrt(max(0, Var_T)),
    Z_Score     = Z_score,
    P_Value     = p_val
  )
  return(list(
    summary = summary_df,
    Sigma = Sigma,
    stats_detail = obs_res$set_stats
  ))
}