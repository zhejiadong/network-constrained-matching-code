#!/usr/bin/env Rscript
# =============================================================================
# run_matching_batch.R
#
# Batch matching for multiple (bz, dep) combinations.
# Files are processed sequentially to avoid CPU contention across Gurobi jobs.
# Within each file, simulations are matched in parallel across cores.
#
# Output: matched_multiprox_bz_{bz}_v{dep}.RData for each combo.
# Skips combos where the output file already exists.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(gurobi)
  library(R.utils)
  library(optmatch)
  library(foreach)
  library(doParallel)
})

# =============================================================================
# CONFIG
# =============================================================================
CODE_ROOT <- "."
DATA_DIR  <- "data"
NET_FILE  <- file.path(DATA_DIR, "net.RData")

source(file.path(CODE_ROOT, "source/utility.R"))
source(file.path(CODE_ROOT, "source/fmnp_gurobi.R"))
source(file.path(CODE_ROOT, "source/run_matching_task.R"))

# Which (bz, dep) combos to run
BZ_VEC  <- 0 # accepted null-only scope
DEP_VEC <- c(1, 2, 3)

# Matching parameters (same as run_matching_parallel_simple.R)
PROX_THRESH  <- 0.5
CALIPER_P1   <- 8
CALIPER_P2   <- 12
TIME_P1      <- 50
TIME_P2      <- 150
GAP_P1       <- 0.01
GAP_P2       <- 0.02
MIN_CTRL     <- 0.01
MAX_CTRL     <- 100
FMNP_THREADS <- 1
N_CORES      <- max(1, parallel::detectCores() - 1)
BATCH_SIZE   <- 20   # print progress every this many sims
FMLA         <- as.formula(Z ~ X1 + X2 + X3)

# =============================================================================
# HELPER: match one simulation
# (identical logic to run_matching_parallel_simple.R)
# =============================================================================
match_one_sim <- function(i, df.list, fmla, net_dist, net_proximity,
                          prox_thresh, time1, time2, gap1, gap2,
                          fmnp_threads) {
  current_df <- df.list[[i]]

  maha.dist_mat <- maha_my(fmla, df = current_df, mode = "tc", cov.type = "pooled")

  fm1 <- tryCatch(
    optmatch::fullmatch(maha.dist_mat, data = current_df,
                        min.controls = MIN_CTRL, max.controls = MAX_CTRL),
    error = function(e) NULL
  )
  current_df$subclass_full1 <- if (!is.null(fm1)) factor(as.numeric(factor(fm1))) else NA

  D_net <- network_penalty(maha.dist_mat, net_dist = net_dist,
                           prox_thresh = prox_thresh, lessequal = FALSE,
                           treat = current_df$Z, penalty = Inf)
  fm2 <- tryCatch(
    optmatch::fullmatch(D_net, data = current_df,
                        min.controls = MIN_CTRL, max.controls = MAX_CTRL),
    error = function(e) NULL
  )
  current_df$subclass_full2 <- if (!is.null(fm2)) factor(as.numeric(factor(fm2))) else NA

  match_res <- run_matching_task(
    data = current_df, fmla = fmla,
    maha.dist_mat = maha.dist_mat,
    net_dist = net_dist, net_proximity = net_proximity,
    prox_thresh = prox_thresh,
    time1 = time1, time2 = time2,
    cal1 = CALIPER_P1, cal2 = CALIPER_P2,
    gap1 = gap1, gap2 = gap2,
    min.controls = MIN_CTRL,
    i_sim = i, threads = fmnp_threads
  )

  if (!is.null(match_res) && !is.null(match_res$data)) {
    current_df$subclass_p_2 <- match_res$data$subclass
    status <- "success"
  } else {
    current_df$subclass_p_2 <- NA
    status <- "failed"
  }

  list(index = i, data = current_df, status = status)
}

# =============================================================================
# MAIN: loop over (bz, dep) sequentially
# =============================================================================
load(NET_FILE)   # net_dist, net_proximity

total_start <- Sys.time()

for (bz in BZ_VEC) {
  for (dep in DEP_VEC) {

    input_file  <- file.path(DATA_DIR, sprintf("bz%g_v%d.RData", bz, dep))
    output_file <- file.path(DATA_DIR, sprintf("matched_multiprox_bz_%g_v%d.RData", bz, dep))

    cat(sprintf("\n========================================\n"))
    cat(sprintf("bz = %g  |  dep = %d\n", bz, dep))
    cat(sprintf("Input : %s\n", input_file))
    cat(sprintf("Output: %s\n", output_file))
    cat(sprintf("========================================\n"))

    if (file.exists(output_file)) {
      cat("Output already exists — skipping.\n")
      next
    }

    if (!file.exists(input_file)) {
      cat("Input file not found — skipping.\n")
      next
    }

    load(input_file)   # df.list (and possibly other objects)
    n_sim <- length(df.list)
    cat(sprintf("Simulations: %d  |  Cores: %d\n", n_sim, N_CORES))

    cl <- makeCluster(N_CORES)
    registerDoParallel(cl)

    clusterExport(cl,
      c("df.list", "net_dist", "net_proximity", "FMLA",
        "PROX_THRESH", "TIME_P1", "TIME_P2", "GAP_P1", "GAP_P2",
        "FMNP_THREADS", "MIN_CTRL", "MAX_CTRL",
        "CALIPER_P1", "CALIPER_P2",
        "match_one_sim", "maha_my", "network_penalty", "run_matching_task",
        "CODE_ROOT"),
      envir = environment()
    )

    clusterEvalQ(cl, {
      library(dplyr); library(Matrix); library(gurobi)
      library(R.utils); library(optmatch)
      source(file.path(CODE_ROOT, "source/fmnp_gurobi.R"))
      NULL
    })

    t0 <- Sys.time()
    make_batches <- function(ids, size) split(ids, ceiling(seq_along(ids) / size))
    batches <- make_batches(seq_len(n_sim), BATCH_SIZE)

    df.final.list <- vector("list", n_sim)
    failure_count <- 0

    for (b in seq_along(batches)) {
      batch_ids <- batches[[b]]
      bt <- Sys.time()

      chunk <- foreach(i = batch_ids,
                       .packages = c("optmatch", "Matrix", "gurobi", "R.utils"),
                       .errorhandling = "pass") %dopar% {
        match_one_sim(i, df.list, FMLA, net_dist, net_proximity,
                      PROX_THRESH, TIME_P1, TIME_P2, GAP_P1, GAP_P2,
                      FMNP_THREADS)
      }

      for (res in chunk) {
        if (!is.list(res) || is.null(res$index)) {
          cat(sprintf("  [ERROR] foreach returned non-list: %s\n",
                      paste(as.character(res), collapse = " ")))
          next
        }
        df.final.list[[res$index]] <- res$data
        if (res$status == "failed") {
          failure_count <- failure_count + 1
          cat(sprintf("  [FAIL] sim %d\n", res$index))
        }
      }

      done <- min(max(batch_ids), n_sim)
      elapsed <- as.numeric(Sys.time() - t0, units = "mins")
      eta     <- if (done > 0) elapsed / done * (n_sim - done) else NA
      cat(sprintf("  [%d/%d] batch %d/%d | +%.1f min elapsed | ETA %.1f min\n",
                  done, n_sim, b, length(batches),
                  as.numeric(Sys.time() - bt, units = "mins"),
                  eta))
    }

    stopCluster(cl)
    cat(sprintf("Done: %d/%d success | total %.1f min\n",
                n_sim - failure_count, n_sim,
                as.numeric(Sys.time() - t0, units = "mins")))

    # Save same objects as run_matching_parallel_simple.R:
    # df.final.list + net objects + everything loaded from the input file
    exclude <- c("cl", "t0", "bt", "b", "batches", "batch_ids",
                 "chunk", "res", "done", "elapsed", "eta",
                 "failure_count", "save_env", "exclude", "obj",
                 "input_file", "output_file", "total_start",
                 "make_batches", "match_one_sim",
                 "BZ_VEC", "DEP_VEC", "PROX_THRESH", "CALIPER_P1",
                 "CALIPER_P2", "TIME_P1", "TIME_P2", "GAP_P1", "GAP_P2",
                 "MIN_CTRL", "MAX_CTRL", "FMNP_THREADS", "N_CORES",
                 "BATCH_SIZE", "FMLA", "CODE_ROOT", "DATA_DIR", "NET_FILE")
    vars_to_save <- setdiff(ls(), exclude)
    save_env <- new.env(parent = emptyenv())
    for (obj in vars_to_save) assign(obj, get(obj), envir = save_env)
    save(list = vars_to_save, file = output_file, envir = save_env)
    cat(sprintf("Saved: %s\n", output_file))
  }
}

cat(sprintf("\n=== All done in %.1f min ===\n",
            as.numeric(Sys.time() - total_start, units = "mins")))
