#!/usr/bin/env Rscript
# =============================================================================
# run_ri_power.R
#
# Power analysis: run RI at the selected (eta*, rho*) for each dep_index,
# across multiple bz values (treatment effects).
#
# For each (dep_index, bz) combination, computes per-sim p-values for:
#   - LM_Crude    : lm(Y ~ Z)  — no matching adjustment
#   - LM_Subclass : lm(Y ~ Z + factor(subclass))  — per matching method
#   - RI_Naive    : RI assuming independent matched sets
#   - RI_Decay    : RI with decay at (eta*, rho*) from ri_optimal_params.csv
#
# Requires:
#   ri_optimal_params.csv  (output of select_optimal_params.R)
#   matched_multiprox_bz_{bz}_v{dep}.RData
#
# Output:
#   ri_power_results.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(foreach)
  library(doParallel)
})

# =============================================================================
# CONFIG
# =============================================================================
CODE_ROOT <- "."
DATA_DIR  <- "data"
GRID_DIR  <- file.path(DATA_DIR, "ri_grid_outputs")
OUT_DIR   <- file.path(DATA_DIR, "power_outputs")

PARAMS_CSV <- file.path(GRID_DIR, "ri_optimal_params.csv")  # from select_optimal_params.R
OUT_CSV    <- file.path(OUT_DIR, "ri_power_results.csv")

# Which (dep_index, bz) combos to run
DEP_VEC <- c(1, 2, 3)     # 1=Independent, 2=Moderate, 3=Strong
BZ_VEC  <- c(0, 0.1, 0.2, 0.3, 0.4)

SUBCLASS_LIST <- list(
  list(name = "Full_1", col = "subclass_full1"),
  list(name = "Full_2", col = "subclass_full2"),
  list(name = "P_2",    col = "subclass_p_2")
)

D0          <- 2
WEIGHT_TYPE <- "ns"
BATCH_SIZE  <- 20
N_CORES     <- max(1, parallel::detectCores() - 1)

# =============================================================================
# SOURCE
# =============================================================================
source(file.path(CODE_ROOT, "source/covariance.R"))
source(file.path(CODE_ROOT, "source/statistics.R"))
source(file.path(CODE_ROOT, "source/utility_ri.R"))

# =============================================================================
# SAFETY CHECK
# =============================================================================
if (file.exists(OUT_CSV)) {
  stop(sprintf(
    "[run_ri_power] Output file already exists:\n  %s\nDelete it before re-running.",
    OUT_CSV))
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# LOAD OPTIMAL PARAMS
# =============================================================================
if (!file.exists(PARAMS_CSV)) {
  stop(sprintf("Optimal params file not found: %s\nRun select_optimal_params.R first.", PARAMS_CSV))
}
optimal <- read.csv(PARAMS_CSV, stringsAsFactors = FALSE)
cat("[run_ri_power] Optimal params:\n")
print(optimal, row.names = FALSE)

# =============================================================================
# HELPER: run one simulation
# =============================================================================
run_one_sim <- function(i, df.final.list, net_dist,
                        SUBCLASS_LIST, D0, WEIGHT_TYPE,
                        eta_star, rho_star, dep_index, bz) {

  fail_row <- function(sc_name, method, msg) data.frame(
    dep_index = dep_index, bz = bz, sim_id = i,
    subclass = sc_name, method = method,
    p_val = NA_real_, status = msg,
    stringsAsFactors = FALSE
  )

  df_sim  <- df.final.list[[i]]
  rows    <- list()

  # LM Crude (no matching) — same for all subclasses
  lm_crude_p <- tryCatch({
    m <- lm(Y ~ Z, data = df_sim)
    summary(m)$coefficients["Z", "Pr(>|t|)"]
  }, error = function(e) NA_real_)

  rows[[length(rows) + 1]] <- data.frame(
    dep_index = dep_index, bz = bz, sim_id = i,
    subclass = "All", method = "LM_Crude",
    p_val = lm_crude_p, status = "ok",
    stringsAsFactors = FALSE
  )

  for (sc in SUBCLASS_LIST) {
    col <- sc$col
    if (!col %in% names(df_sim) || all(is.na(df_sim[[col]]))) {
      rows[[length(rows) + 1]] <- fail_row(sc$name, "all", "missing_subclass")
      next
    }
    df_sub          <- df_sim[!is.na(df_sim[[col]]), ]
    df_sub$subclass <- df_sub[[col]]

    # --- LM Subclass ---
    lm_sub_p <- tryCatch({
      m <- lm(Y ~ Z + factor(subclass), data = df_sub)
      summary(m)$coefficients["Z", "Pr(>|t|)"]
    }, error = function(e) NA_real_)

    rows[[length(rows) + 1]] <- data.frame(
      dep_index = dep_index, bz = bz, sim_id = i,
      subclass = sc$name, method = "LM_Subclass",
      p_val = lm_sub_p, status = "ok",
      stringsAsFactors = FALSE
    )

    # --- RI Naive ---
    ri_naive <- tryCatch(
      run_randomization_inference(df_sub, cov_type = "naive",
                                  weight_type = WEIGHT_TYPE,
                                  subclass_col = "subclass"),
      error = function(e) NULL
    )
    rows[[length(rows) + 1]] <- if (!is.null(ri_naive)) {
      data.frame(dep_index = dep_index, bz = bz, sim_id = i,
                 subclass = sc$name, method = "RI_Naive",
                 p_val = ri_naive$summary$P_Value, status = "ok",
                 stringsAsFactors = FALSE)
    } else {
      fail_row(sc$name, "RI_Naive", "ERROR:naive")
    }

    # --- RI Decay at (eta*, rho*) — precompute for speed ---
    ri_decay_p <- tryCatch({
      obs_res    <- get_observed_stats(df_sub, WEIGHT_TYPE, "subclass", "Z", "Y")
      design_cov <- precompute_design_cov_matrix(df_sub, "subclass", "Z")
      set_dist_m <- get_set_distance_matrix(net_dist, df_sub$subclass, "dmin")
      wmat       <- get_set_decay_matrix(set_dist_m, eta_star, rho_star, D0, 1)
      Sigma      <- get_U_cov_from_precomputed(design_cov, wmat)
      wvec       <- obs_res$set_stats$weight
      Var_T      <- as.numeric(t(wvec) %*% Sigma %*% wvec)
      if (!is.finite(Var_T) || Var_T <= 0) stop("non-positive variance")
      Z_sc       <- (obs_res$T_obs - obs_res$E_T) / sqrt(Var_T)
      2 * (1 - pnorm(abs(Z_sc)))
    }, error = function(e) NA_real_)

    rows[[length(rows) + 1]] <- data.frame(
      dep_index = dep_index, bz = bz, sim_id = i,
      subclass = sc$name, method = "RI_Decay",
      p_val = ri_decay_p, status = if (is.na(ri_decay_p)) "ERROR:decay" else "ok",
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

# =============================================================================
# MAIN LOOP over (dep_index, bz)
# =============================================================================
fns_to_export <- c(
  "get_design_cov_bound", "get_U_covariance_matrix",
  "get_set_distance_matrix", "get_set_decay_matrix",
  "get_set_mean_decay_matrix", "get_set_mean_decay_matrix_soft",
  "precompute_design_cov_matrix", "get_U_cov_from_precomputed",
  "get_observed_stats", "get_U_moments",
  "run_randomization_inference"
)

make_batches <- function(ids, size) split(ids, ceiling(seq_along(ids) / size))

for (dep in DEP_VEC) {
  for (bz in BZ_VEC) {

    matched_file <- file.path(DATA_DIR,
      sprintf("matched_multiprox_bz_%g_v%d.RData", bz, dep))

    if (!file.exists(matched_file)) {
      cat(sprintf("[run_ri_power] SKIP dep=%d bz=%g — file not found\n", dep, bz))
      next
    }

    # Get optimal params for this dep
    opt_row <- optimal[optimal$dep_index == dep, ]
    if (nrow(opt_row) == 0) {
      cat(sprintf("[run_ri_power] SKIP dep=%d — no optimal params\n", dep))
      next
    }
    eta_star <- opt_row$eta[1]
    rho_star <- opt_row$rho[1]

    cat(sprintf("\n[run_ri_power] dep=%d  bz=%g  eta*=%.2f  rho*=%.2f\n",
                dep, bz, eta_star, rho_star))

    load(matched_file)   # df.final.list
    load(file.path(DATA_DIR, "net.RData"))   # net_dist

    n_rep   <- length(df.final.list)
    batches <- make_batches(seq_len(n_rep), BATCH_SIZE)

    cl <- makeCluster(N_CORES)
    registerDoParallel(cl)
    clusterExport(cl,
      c("df.final.list", "net_dist", "SUBCLASS_LIST", "D0", "WEIGHT_TYPE",
        "eta_star", "rho_star", "dep", "bz",
        "run_one_sim", fns_to_export),
      envir = environment()
    )

    for (b in seq_along(batches)) {
      batch_ids <- batches[[b]]
      cat(sprintf("  Batch %d/%d (sims %d-%d)\n",
                  b, length(batches), min(batch_ids), max(batch_ids)))

      chunk <- foreach(i = batch_ids, .combine = rbind,
                       .packages = c("stats", "dplyr")) %dopar% {
        tryCatch(
          run_one_sim(i, df.final.list, net_dist,
                      SUBCLASS_LIST, D0, WEIGHT_TYPE,
                      eta_star, rho_star, dep, bz),
          error = function(e) data.frame(
            dep_index = dep, bz = bz, sim_id = i,
            subclass = "all", method = "all",
            p_val = NA_real_,
            status = paste("ERROR:", conditionMessage(e)),
            stringsAsFactors = FALSE
          )
        )
      }

      ok <- chunk[!is.na(chunk$status) & chunk$status == "ok", , drop = FALSE]
      if (nrow(ok) > 0) {
        write.table(ok[, setdiff(names(ok), "status")],
                    OUT_CSV, sep = ",", row.names = FALSE,
                    col.names = !file.exists(OUT_CSV),
                    append    = file.exists(OUT_CSV),
                    quote     = FALSE)
      }
    }

    stopCluster(cl)
  }
}

# =============================================================================
# SUMMARY
# =============================================================================
results <- read.csv(OUT_CSV, stringsAsFactors = FALSE)
cat("\n--- Rejection rates (alpha = 0.05) ---\n")
summary_tbl <- results %>%
  group_by(dep_index, bz, subclass, method) %>%
  summarise(reject_rate = mean(p_val < 0.05, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  arrange(dep_index, bz, subclass, method)
print(as.data.frame(summary_tbl), digits = 3, row.names = FALSE)

cat(sprintf("\n[run_ri_power] Results -> %s\n", OUT_CSV))
cat("[run_ri_power] Done.\n")
