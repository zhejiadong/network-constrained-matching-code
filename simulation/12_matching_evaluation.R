#!/usr/bin/env Rscript
# =============================================================================
# run_matching_eval.R
#
# Compute three matching evaluation metrics across all simulation reps at bz=0
# for dep_index = 1 (Independent), 2 (Moderate), 3 (Strong).
#
# Metrics:
#   (1) Covariate balance  — |SMD| before vs after matching
#   (2) Covariate similarity — average within-set Mahalanobis distance
#   (3) Network separation   — average within-set pairwise network distance
#
# Outputs (-> diagnosis/):
#   matching_eval_dep3.pdf   main manuscript: 3-panel figure (dep=3 only)
#   matching_eval_all.pdf    appendix: same 3 panels faceted by dep_index
#
# Sources: evaluation/matching_metrics.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(foreach)
  library(doParallel)
})

# =============================================================================
# CONFIG
# =============================================================================
CODE_ROOT <- "."
DATA_DIR  <- "data"
OUT_DIR   <- file.path(DATA_DIR, "diagnosis")

DEP_VEC <- c(1, 2, 3)
BZ      <- 0

COVARS <- c("X1", "X2", "X3")

SUBCLASS_METHODS <- list(
  "Full_1" = "subclass_full1",
  "Full_2" = "subclass_full2",
  "P_2"    = "subclass_p_2"
)

DEP_LABELS   <- c("1" = "Independent", "2" = "Moderate", "3" = "Strong")
METHOD_ORDER <- c("Before", "Full_1", "Full_2", "P_2")
METHOD_LABELS <- c(
  "Before" = "Before matching",
  "Full_1" = "Covariates only",
  "Full_2" = "Single penalty",
  "P_2"    = "Dual penalty"
)
METHOD_COLORS <- c(
  "Before" = "#F8766D",   # red (matches reference)
  "Full_1" = "#7CAE00",   # olive green
  "Full_2" = "#00BFC4",   # teal/cyan
  "P_2"    = "#C77CFF"    # purple
)

N_CORES <- max(1, parallel::detectCores() - 1)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
source(file.path(CODE_ROOT, "source/matching_metrics.R"))

# Functions to export to parallel workers
FNS_EXPORT <- c(
  "compute_smd_long",  "smd_before",    "smd_after",
  "compute_mahal_long","mahal_before",  "mahal_after",  ".get_S_inv",
  "compute_net_long",  "net_before_stats", "net_pair_stats", "net_across_dmin_stats"
)

# =============================================================================
# COLLECT METRICS ACROSS SIMS AND DEP_INDEX  (skip if cache exists)
# =============================================================================
cache_files <- file.path(OUT_DIR, c("eval_smd_agg.rds", "eval_df_mahal.rds",
                                    "eval_df_net_plot.rds", "eval_comp_agg.rds",
                                    "eval_dmin_agg.rds"))

if (all(file.exists(cache_files))) {
  cat("[run_matching_eval] Cache found — loading prepared tables, skipping computation.\n")
  smd_agg     <- readRDS(cache_files[1])
  df_mahal    <- readRDS(cache_files[2])
  df_net_plot <- readRDS(cache_files[3])
  comp_agg    <- readRDS(cache_files[4])
  dmin_agg    <- readRDS(cache_files[5])
} else {

smd_all   <- list()
mahal_all <- list()
net_all   <- list()

for (dep in DEP_VEC) {

  matched_file <- file.path(DATA_DIR,
    sprintf("matched_multiprox_bz_%g_v%d.RData", BZ, dep))
  net_file <- file.path(DATA_DIR, "net.RData")

  if (!file.exists(matched_file)) {
    cat(sprintf("[run_matching_eval] SKIP dep=%d — file not found\n", dep))
    next
  }

  cat(sprintf("[run_matching_eval] dep=%d: loading data...\n", dep))
  load(matched_file)   # df.final.list
  load(net_file)       # net_dist

  n_rep <- length(df.final.list)
  cat(sprintf("  %d sims | network %dx%d | %d cores\n",
              n_rep, nrow(net_dist), ncol(net_dist), N_CORES))

  cl <- makeCluster(N_CORES)
  registerDoParallel(cl)
  clusterExport(cl,
    c("df.final.list", "net_dist", "COVARS", "SUBCLASS_METHODS", "dep",
      FNS_EXPORT),
    envir = environment()
  )

  results <- foreach(i = seq_len(n_rep),
                     .packages = c("dplyr", "MASS"),
                     .errorhandling = "pass") %dopar% {
    df_sim <- df.final.list[[i]]
    list(
      smd   = cbind(dep_index = dep, sim_id = i,
                    compute_smd_long(df_sim, COVARS, SUBCLASS_METHODS),
                    stringsAsFactors = FALSE),
      mahal = cbind(dep_index = dep, sim_id = i,
                    compute_mahal_long(df_sim, COVARS, SUBCLASS_METHODS),
                    stringsAsFactors = FALSE),
      net   = cbind(dep_index = dep, sim_id = i,
                    compute_net_long(df_sim, net_dist, SUBCLASS_METHODS),
                    stringsAsFactors = FALSE)
    )
  }
  stopCluster(cl)

  # Filter out any error returns
  ok <- Filter(is.list, results)
  cat(sprintf("  %d/%d sims completed successfully\n", length(ok), n_rep))

  smd_all[[dep]]   <- bind_rows(lapply(ok, `[[`, "smd"))
  mahal_all[[dep]] <- bind_rows(lapply(ok, `[[`, "mahal"))
  net_all[[dep]]   <- bind_rows(lapply(ok, `[[`, "net"))
}

df_smd   <- bind_rows(smd_all)
df_mahal <- bind_rows(mahal_all)
df_net   <- bind_rows(net_all)

# =============================================================================
# AGGREGATE
# =============================================================================

# SMD: mean + IQR across sims per (dep_index, method, covariate)
smd_agg <- df_smd %>%
  group_by(dep_index, method, covariate) %>%
  summarise(
    mean_smd = mean(abs_smd,             na.rm = TRUE),
    q25_smd  = quantile(abs_smd, 0.25,   na.rm = TRUE),
    q75_smd  = quantile(abs_smd, 0.75,   na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(
    dep_label = factor(DEP_LABELS[as.character(dep_index)],
                       levels = c("Independent", "Moderate", "Strong")),
    method    = factor(method, levels = METHOD_ORDER)
  )

# Mahal: keep per-sim values for boxplot
df_mahal <- df_mahal %>%
  mutate(
    dep_label = factor(DEP_LABELS[as.character(dep_index)],
                       levels = c("Independent", "Moderate", "Strong")),
    method    = factor(method, levels = METHOD_ORDER)
  )

# Net: within-set (and "all" for Before) for boxplot
df_net_plot <- df_net %>%
  filter(scope %in% c("all", "within")) %>%
  mutate(
    dep_label = factor(DEP_LABELS[as.character(dep_index)],
                       levels = c("Independent", "Moderate", "Strong")),
    method    = factor(method, levels = METHOD_ORDER)
  )

# Net composition: mean proportions of d=1 / d=2 / d>=3 across sims,
# within uses all pairwise within-set; across uses all pairwise across-set.
comp_agg <- df_net %>%
  filter(scope %in% c("all", "within", "across")) %>%
  group_by(dep_index, method, scope) %>%
  summarise(
    mean_p1   = mean(p1,   na.rm = TRUE),
    mean_p2   = mean(p2,   na.rm = TRUE),
    mean_pge3 = mean(pge3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols = c(mean_p1, mean_p2, mean_pge3),
                      names_to  = "band",
                      values_to = "proportion") %>%
  mutate(
    band = factor(
      recode(band, mean_p1 = "d = 1", mean_p2 = "d = 2", mean_pge3 = "d >= 3"),
      levels = c("d >= 3", "d = 2", "d = 1")   # stack order: d=1 on bottom
    ),
    dep_label = factor(DEP_LABELS[as.character(dep_index)],
                       levels = c("Independent", "Moderate", "Strong")),
    method    = factor(method, levels = METHOD_ORDER)
  )

# dmin composition: across-set set-pair minimum distances
dmin_agg <- df_net %>%
  filter(scope == "across_dmin", method != "Before") %>%
  group_by(dep_index, method) %>%
  summarise(
    mean_p1   = mean(p1,   na.rm = TRUE),
    mean_p2   = mean(p2,   na.rm = TRUE),
    mean_pge3 = mean(pge3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols = c(mean_p1, mean_p2, mean_pge3),
                      names_to  = "band",
                      values_to = "proportion") %>%
  mutate(
    band = factor(
      recode(band, mean_p1 = "d = 1", mean_p2 = "d = 2", mean_pge3 = "d >= 3"),
      levels = c("d >= 3", "d = 2", "d = 1")
    ),
    dep_label = factor(DEP_LABELS[as.character(dep_index)],
                       levels = c("Independent", "Moderate", "Strong")),
    method    = factor(method, levels = METHOD_ORDER)
  )

# =============================================================================
# CACHE PREPARED TABLES (consumed by slide_eval_cov_comp.R and similar)
# =============================================================================
saveRDS(smd_agg,     file.path(OUT_DIR, "eval_smd_agg.rds"))
saveRDS(df_mahal,    file.path(OUT_DIR, "eval_df_mahal.rds"))
saveRDS(df_net_plot, file.path(OUT_DIR, "eval_df_net_plot.rds"))
saveRDS(comp_agg,    file.path(OUT_DIR, "eval_comp_agg.rds"))
saveRDS(dmin_agg,    file.path(OUT_DIR, "eval_dmin_agg.rds"))
cat("[run_matching_eval] Cached prepared tables -> eval_*.rds\n")

} # end else (compute from scratch)

# =============================================================================
# PLOT BUILDERS
# =============================================================================

make_love_plot <- function(data, show_title = TRUE) {
  dodge <- position_dodge(width = 0.5)
  p <- ggplot(data,
              aes(x = mean_smd, y = covariate,
                  colour = method, group = method)) +
    geom_vline(xintercept = 0.1, linetype = "dashed",
               colour = "grey50", linewidth = 0.6) +
    geom_errorbarh(aes(xmin = q25_smd, xmax = q75_smd),
                   height = 0.15, linewidth = 0.8, alpha = 0.7,
                   position = dodge) +
    geom_point(size = 3, position = dodge) +
    scale_colour_manual(values = METHOD_COLORS,
                        labels = METHOD_LABELS,
                        name   = "Method") +
    scale_x_continuous(name = "Absolute SMD") +
    scale_y_discrete(name = "Covariate") +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          legend.position  = "right")
  if (show_title)
    p <- p + labs(title = "Covariate Balance (SMD)")
  p
}

make_mahal_plot <- function(data, show_title = TRUE) {
  p <- ggplot(data, aes(x = method, y = mahal_dist, fill = method)) +
    geom_boxplot(outlier.size = 0.7, outlier.alpha = 0.4, width = 0.55) +
    scale_fill_manual(values = METHOD_COLORS, guide = "none") +
    scale_x_discrete(labels = METHOD_LABELS, name = "Matching method") +
    scale_y_continuous(name = "Average distance") +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          axis.text.x      = element_text(size = 11))
  if (show_title)
    p <- p + labs(title = expression("Covariate Similarity (" * delta[ij] * ")"))
  p
}

make_net_plot <- function(data, show_title = TRUE) {
  p <- ggplot(data, aes(x = method, y = avg, fill = method)) +
    geom_boxplot(outlier.size = 0.7, outlier.alpha = 0.4, width = 0.55) +
    scale_fill_manual(values = METHOD_COLORS, guide = "none") +
    scale_x_discrete(labels = METHOD_LABELS, name = "Matching method") +
    scale_y_continuous(name = "Average distance") +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          axis.text.x      = element_text(size = 11))
  if (show_title)
    p <- p + labs(title = expression("Network Proximity (" * d[ij] * ")"))
  p
}

BAND_COLORS <- c("d = 1" = "#D73027", "d = 2" = "#FC8D59", "d >= 3" = "#FEE08B")

make_comp_plot <- function(data, show_title = TRUE,
                           scopes = c("within", "across")) {
  plot_df <- data %>%
    filter(method == "Before" | scope %in% scopes) %>%
    mutate(
      method_scope = dplyr::case_when(
        method == "Before" ~ "Before",
        scope == "within"  ~ paste0(as.character(method), "_within"),
        scope == "across"  ~ paste0(as.character(method), "_across"),
        TRUE               ~ as.character(method)
      )
    )

  canonical_levels <- c("Before",
                        "Full_1_within", "Full_1_across",
                        "Full_2_within", "Full_2_across",
                        "P_2_within",    "P_2_across")
  levels_keep <- canonical_levels[canonical_levels %in% unique(plot_df$method_scope)]
  plot_df$method_scope <- factor(plot_df$method_scope, levels = levels_keep)

  # If only one of within/across is shown, drop the suffix from the labels.
  suffix_needed <- length(scopes) > 1
  method_scope_labels <- c(
    "Before"        = "Before matching",
    "Full_1_within" = if (suffix_needed) "Cov. only (within)" else "Covariates only",
    "Full_1_across" = if (suffix_needed) "Cov. only (across)" else "Covariates only",
    "Full_2_within" = if (suffix_needed) "Single pen. (within)" else "Single penalty",
    "Full_2_across" = if (suffix_needed) "Single pen. (across)" else "Single penalty",
    "P_2_within"    = if (suffix_needed) "Dual pen. (within)"   else "Dual penalty",
    "P_2_across"    = if (suffix_needed) "Dual pen. (across)"   else "Dual penalty"
  )

  p <- ggplot(plot_df, aes(x = method_scope, y = proportion, fill = band)) +
    geom_col(width = 0.6, colour = "grey30", linewidth = 0.25) +
    scale_fill_manual(values = BAND_COLORS,
                      name   = "Distance",
                      guide  = guide_legend(reverse = TRUE)) +
    scale_x_discrete(labels = method_scope_labels, name = "Matching method") +
    scale_y_continuous(labels = scales::percent, name = "Proportion") +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          axis.text.x      = element_text(size = 11))
  if (show_title)
    p <- p + labs(title = "Distance Composition: Before and after matching")
  p
}

make_dmin_plot <- function(data, show_title = TRUE) {
  p <- ggplot(data, aes(x = method, y = proportion, fill = band)) +
    geom_col(width = 0.6, colour = "grey30", linewidth = 0.25) +
    scale_fill_manual(values = BAND_COLORS,
                      name   = "Distance",
                      guide  = guide_legend(reverse = TRUE)) +
    scale_x_discrete(labels = METHOD_LABELS, name = "Matching method") +
    scale_y_continuous(labels = scales::percent, name = "Proportion") +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          axis.text.x      = element_text(size = 11))
  if (show_title)
    p <- p + labs(title = "Distance Composition: across-set minimum distances")
  p
}

# =============================================================================
# OUTPUT 1: dep=3 only — main manuscript
# =============================================================================
cat("\n[run_matching_eval] Building main manuscript figures (dep=3)...\n")

p1 <- make_love_plot( smd_agg     %>% filter(dep_index == 3))
p2 <- make_mahal_plot(df_mahal    %>% filter(dep_index == 3))
p3 <- make_net_plot(  df_net_plot %>% filter(dep_index == 3))
p4 <- make_comp_plot( comp_agg    %>% filter(dep_index == 3))
p5 <- make_dmin_plot( dmin_agg    %>% filter(dep_index == 3))

# Figure A: covariate metrics (2 panels)
p_cov <- p1 | p2
ggsave(file.path(OUT_DIR, "matching_eval_dep3_cov.pdf"),
       plot = p_cov, width = 12, height = 5)
cat("[run_matching_eval] Saved -> matching_eval_dep3_cov.pdf\n")

# Figure B: network metrics (2 panels)
p_net <- p3 | p4
ggsave(file.path(OUT_DIR, "matching_eval_dep3_net.pdf"),
       plot = p_net, width = 15, height = 6)
cat("[run_matching_eval] Saved -> matching_eval_dep3_net.pdf\n")

# Figure B': covariate similarity | network proximity (main, dep=3)
p_simnet <- p2 | p3
ggsave(file.path(OUT_DIR, "matching_eval_dep3_simnet.pdf"),
       plot = p_simnet, width = 12, height = 5)
cat("[run_matching_eval] Saved -> matching_eval_dep3_simnet.pdf\n")

# Figure B'': network proximity | within-only composition (main, dep=3)
p4_within <- make_comp_plot(comp_agg %>% filter(dep_index == 3),
                            scopes = "within")
p_netcomp <- p3 | p4_within
ggsave(file.path(OUT_DIR, "matching_eval_dep3_netcomp.pdf"),
       plot = p_netcomp, width = 13, height = 5)
cat("[run_matching_eval] Saved -> matching_eval_dep3_netcomp.pdf\n")

# Figure C: dmin (appendix standalone, dep=3)
ggsave(file.path(OUT_DIR, "diag_across_dmin_barplot.pdf"),
       plot = p5, width = 7, height = 6)
cat("[run_matching_eval] Saved -> diag_across_dmin_barplot.pdf\n")

# =============================================================================
# OUTPUT 2: all dep_index — appendix (faceted rows, two separate PDFs)
# =============================================================================
cat("[run_matching_eval] Building appendix figures (all dep)...\n")

p1_all <- make_love_plot( smd_agg,     show_title = FALSE) +
            facet_wrap(~ dep_label, ncol = 3) +
            labs(title = "Covariate Balance: |SMD|")

TWO_LINE_LABELS <- c("Before" = "Before\nmatching",
                     "Full_1" = "Covariates\nonly",
                     "Full_2" = "Single\npenalty",
                     "P_2"    = "Dual\npenalty")

p2_all <- make_mahal_plot(df_mahal,    show_title = FALSE) +
            facet_wrap(~ dep_label, ncol = 3) +
            labs(title = expression("Covariate Similarity (" * delta[ij] * ")")) +
            scale_x_discrete(labels = TWO_LINE_LABELS, name = "Matching method")

p3_all <- make_net_plot(  df_net_plot, show_title = FALSE) +
            facet_wrap(~ dep_label, ncol = 3) +
            labs(title = expression("Network Proximity (" * d[ij] * ")")) +
            scale_x_discrete(labels = TWO_LINE_LABELS, name = "Matching method")

p4_all <- make_comp_plot( comp_agg,    show_title = FALSE, scopes = "within") +
            facet_wrap(~ dep_label, ncol = 3) +
            labs(title = expression("Within-set " * d[ij] * " distribution")) +
            scale_x_discrete(labels = c("Before"        = "Before\nmatching",
                                        "Full_1_within" = "Covariates\nonly",
                                        "Full_2_within" = "Single\npenalty",
                                        "P_2_within"    = "Dual\npenalty"),
                             name = "Matching method")

p5_all <- make_dmin_plot( dmin_agg,   show_title = FALSE) +
            facet_wrap(~ dep_label, ncol = 3) +
            labs(title = "Across-Set Minimum Distance Distribution")

# Combined 3-row figure: covariate similarity / network proximity / distance composition
p_appendix_all3 <- p2_all / p3_all / p4_all

ggsave(file.path(OUT_DIR, "matching_eval_all_combined.pdf"),
       plot = p_appendix_all3, width = 15, height = 15)
cat("[run_matching_eval] Saved -> matching_eval_all_combined.pdf\n")

# Appendix B': covariate similarity (top) over network proximity (bottom)
p_appendix_simnet <- p2_all / p3_all +
  plot_annotation(title = "Covariate Similarity & Network Proximity — All Dependence Levels (bz = 0)")

ggsave(file.path(OUT_DIR, "matching_eval_all_simnet.pdf"),
       plot = p_appendix_simnet, width = 15, height = 10)
cat("[run_matching_eval] Saved -> matching_eval_all_simnet.pdf\n")

# Appendix C: dmin across all dep levels
ggsave(file.path(OUT_DIR, "diag_across_dmin_barplot_all.pdf"),
       plot = p5_all, width = 15, height = 5)
cat("[run_matching_eval] Saved -> diag_across_dmin_barplot_all.pdf\n")

cat("[run_matching_eval] Done.\n")


