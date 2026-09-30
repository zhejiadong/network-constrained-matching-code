#!/usr/bin/env Rscript
# =============================================================================
# select_optimal_params.R
#
# For each dep_index (1=Independent, 2=Moderate, 3=Strong):
#   - Read ri_grid_results_bz0_v{dep}.csv
#   - Filter P_2, method = RI_Decay_dmin
#   - Compute rejection rate per (eta, rho)
#   - Select (eta*, rho*) closest to alpha = 0.05
#
# Outputs:
#   ri_optimal_params.csv   — one row per dep_index: dep, eta*, rho*, reject_rate
#   ri_table2_ri_decay.csv  — rejection rates at (eta*, rho*) for all 3 subclasses
#                             (fills RI Decay rows of Table 2)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

DATA_DIR <- "data"
OUT_DIR  <- file.path(DATA_DIR, "ri_grid_outputs")
ALPHA    <- 0.05
DEP_LABELS <- c("1" = "Independent", "2" = "Moderate", "3" = "Strong")

# ---------------------------------------------------------------------------
# Load all available bz=0 grid CSVs
# ---------------------------------------------------------------------------
csv_files <- list(
  "1" = file.path(OUT_DIR, "ri_grid_results_bz0_v1.csv"),
  "2" = file.path(OUT_DIR, "ri_grid_results_bz0_v2.csv"),
  "3" = file.path(OUT_DIR, "ri_grid_results.csv")          # legacy v3 file
)

df_all <- lapply(names(csv_files), function(dep) {
  f <- csv_files[[dep]]
  if (!file.exists(f)) {
    message(sprintf("  [skip] dep=%s file not found: %s", dep, f))
    return(NULL)
  }
  d <- read.csv(f, stringsAsFactors = FALSE)
  # Back-fill dep_index for old v3 CSV that may lack the column
  if (!"dep_index" %in% names(d)) d$dep_index <- as.integer(dep)
  if (!"bz"        %in% names(d)) d$bz        <- 0
  d
}) %>% Filter(Negate(is.null), .)

if (length(df_all) == 0) stop("No grid CSVs found. Run run_ri_decay_null.R first.")
df_all <- bind_rows(df_all)

cat(sprintf("Loaded %d rows across %d dep_index value(s).\n",
            nrow(df_all), length(unique(df_all$dep_index))))

# ---------------------------------------------------------------------------
# Step 1: Select (eta*, rho*) per dep_index using P_2 rejection rate
# ---------------------------------------------------------------------------
p2_summary <- df_all %>%
  filter(subclass == "P_2", method == "RI_Decay_dmin") %>%
  group_by(dep_index, eta, rho) %>%
  summarise(
    n_sims      = n(),
    reject_rate = mean(p_val < ALPHA, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(dist_to_alpha = abs(reject_rate - ALPHA))

# Only grid-select for Moderate (2) and Strong (3)
# Independent (1): use eta=0, rho=0 (RI_Naive, no decay correction needed)
optimal_grid <- p2_summary %>%
  filter(dep_index %in% c(2, 3)) %>%
  group_by(dep_index) %>%
  slice_min(dist_to_alpha, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(dep_label = DEP_LABELS[as.character(dep_index)]) %>%
  select(dep_index, dep_label, eta, rho, reject_rate, n_sims)

indep_row <- data.frame(
  dep_index   = 1L,
  dep_label   = "Independent",
  eta         = 0,
  rho         = 0,
  reject_rate = NA_real_,
  n_sims      = NA_integer_,
  stringsAsFactors = FALSE
)

optimal <- bind_rows(indep_row, optimal_grid) %>%
  arrange(dep_index)

cat("\n--- Optimal (eta*, rho*) per dep_index (P_2, closest to alpha=0.05) ---\n")
print(as.data.frame(optimal), digits = 3, row.names = FALSE)

out_params <- file.path(OUT_DIR, "ri_optimal_params.csv")
write.csv(optimal, out_params, row.names = FALSE)
cat(sprintf("\nSaved -> %s\n", out_params))

# ---------------------------------------------------------------------------
# Step 2: Table 2 RI Decay fill — rejection rates at (eta*, rho*) for all
#         3 subclasses and all dep_index values
# ---------------------------------------------------------------------------
table2_rows <- lapply(seq_len(nrow(optimal)), function(k) {
  dep  <- optimal$dep_index[k]
  eta_k <- optimal$eta[k]
  rho_k <- optimal$rho[k]

  df_all %>%
    filter(dep_index == dep,
           method == "RI_Decay_dmin",
           abs(eta - eta_k) < 1e-9,
           abs(rho - rho_k) < 1e-9) %>%
    group_by(dep_index, subclass) %>%
    summarise(
      eta         = eta_k,
      rho         = rho_k,
      n_sims      = n(),
      reject_rate = mean(p_val < ALPHA, na.rm = TRUE),
      .groups = "drop"
    )
}) %>% bind_rows()

table2_rows <- table2_rows %>%
  mutate(dep_label = DEP_LABELS[as.character(dep_index)]) %>%
  arrange(dep_index, subclass)

cat("\n--- Table 2 RI Decay fill (at optimal params) ---\n")
print(as.data.frame(table2_rows), digits = 3, row.names = FALSE)

out_table2 <- file.path(OUT_DIR, "ri_table2_ri_decay.csv")
write.csv(table2_rows, out_table2, row.names = FALSE)
cat(sprintf("Saved -> %s\n", out_table2))

cat("\n[select_optimal_params] Done.\n")
