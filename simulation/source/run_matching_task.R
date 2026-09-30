# Main matching task function (two-pass strategy)
# Dependencies: R.utils, fmnp_gurobi.R, utility.R

#' Run matching task with two-pass strategy (strict then relaxed)
#'
#' @param data Data frame for a single simulation.
#' @param fmla Formula for covariates.
#' @param maha.dist_mat Mahalanobis distance matrix (treated × control).
#' @param net_dist Network distance matrix (full N × N).
#' @param net_proximity Network proximity matrix (1 / net_dist).
#' @param prox_thresh Proximity threshold for network penalty.
#' @param time1 Time limit for first pass (seconds).
#' @param time2 Time limit for second pass (seconds).
#' @param gap1 MIP gap for first pass.
#' @param gap2 MIP gap for second pass.
#' @param i_sim Simulation index (for logging).
#' @param threads Number of Gurobi threads to allocate for each FMNP solve.
#'
#' @return The matching result (list with $data and $res) or NULL if both passes fail.
run_matching_task <- function(data, fmla, maha.dist_mat, net_dist, net_proximity,
                              prox_thresh, time1 = 12, time2 = 30, cal1 = 5, cal2 = 8.0,
                              gap1 = 0.01, gap2 = 0.015, min.controls = 0.01, i_sim,
                              threads = 1L) {
  # 1. Define D (Network Penalty) using the dynamic prox_thresh
  D <- network_penalty(maha.dist_mat, 
                       net_dist = net_dist, 
                       prox_thresh = prox_thresh, 
                       lessequal = FALSE, 
                       treat = data$Z, 
                       penalty = Inf)
  # 2. Define parameters for Pass 1 (Strict) and Pass 2 (Relaxed)
  #    Structure: list(caliper, max.controls, timelimit, mipgap, desc)
  passes <- list(
    list(cal = cal1, max_c = 100, time = time1,  gap = gap1,  desc = "Pass 1"),
    list(cal = cal2, max_c = 100,  time = time2 , gap = gap2, desc = "Retry")
  )
  
  final_res <- NULL
  
  for (p in passes) {
    # If we already succeeded in a previous pass, skip
    if (!is.null(final_res)) break 
    
    final_res <- tryCatch({
      R.utils::withTimeout({
        FMNP_gurobi(
          fmla,
          data = data,
          D = D,
          net_proximity = net_proximity,
          caliper = p$cal,
          prox_thresh = prox_thresh,
          min.controls = min.controls,
          max.controls = p$max_c,
          timelimit = p$time,
          mipgap = p$gap,
          thread = threads,
          OutputFlag = 0
        )
      }, timeout = p$time + 15, onTimeout = "error")
    }, error = function(e) {
      # Optional: Only print if it's the Retry pass or if you want detailed logs
      # cat(sprintf("   - %s failed for prox %.2f: %s\n", p$desc, target_prox_thresh, conditionMessage(e)))
      return(NULL)
    })
    
    # Check success
    if (!is.null(final_res)) {
      # You might want to cat() here only if it was a Retry to keep logs clean
      if(p$desc == "Retry") cat(sprintf("   ✅ Sim %d recovered (prox %.2f)\n", i_sim, prox_thresh))
    }
  }
  
  return(final_res)
}