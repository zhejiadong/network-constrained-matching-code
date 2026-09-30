# Export saved accepted P_2 grid rows to the legacy selector's input columns.
# No inference is recomputed; other designs and unavailable statistics are not invented.
files <- file.path("tasks", sprintf("task_%03d.rds", 1:150))
stopifnot(all(file.exists(files)))
objects <- lapply(files, readRDS)
for (i in seq_along(objects)) {
  x <- objects[[i]]
  expected <- (((i - 1L) %% 50L) * 10L + 1L):(((i - 1L) %% 50L) * 10L + 10L)
  stopifnot(x$identity$task == i, x$identity$mc == 500L,
            x$identity$dep == (i - 1L) %/% 50L + 1L,
            nrow(x$grid) == 400L,
            identical(as.integer(x$identity$reps), expected),
            identical(sort(unique(x$grid$sim_id)), expected))
}
grid <- do.call(rbind, lapply(objects, `[[`, "grid"))
stopifnot(nrow(grid) == 60000L, all(grid$design == "P_2"),
          !anyDuplicated(grid[c("dep_index", "sim_id", "eta", "rho")]),
          all(is.finite(grid$p_value)), all(grid$p_value >= 0 & grid$p_value <= 1))
legacy <- data.frame(dep_index = grid$dep_index, bz = 0,
                     sim_id = grid$sim_id, subclass = grid$design,
                     method = "RI_Decay_dmin", eta = grid$eta,
                     rho = grid$rho, p_val = grid$p_value)
out <- file.path("data", "ri_grid_outputs")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
names_out <- c("ri_grid_results_bz0_v1.csv", "ri_grid_results_bz0_v2.csv",
               "ri_grid_results.csv")
for (dep in 1:3) {
  write.csv(legacy[legacy$dep_index == dep, ], file.path(out, names_out[dep]),
            row.names = FALSE)
}
cat("Exported 60000 saved P_2 grid rows; no inference run.\n")
