# Replot-only manuscript assets from unchanged saved simulation outputs.
# Adapted from the existing slide_eval_cov_comp.R, plot_ri_grid.R and plot_power.R.
# Run from simulation/; the optional first argument is a local figure directory.
# Does not run simulations or modify cached data or the original manuscript.
args <- commandArgs(trailingOnly=TRUE)
PLOTS <- if (length(args)) args[1] else file.path('data', 'manuscript_figures')
dir.create(PLOTS, recursive=TRUE, showWarnings=FALSE)
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(patchwork)})

local({
DATA_DIR <- "data"
OUT_DIR  <- file.path(DATA_DIR, "diagnosis")

DEP <- 3   # slide focuses on Strong dependence

METHOD_ORDER  <- c("Before", "Full_1", "Full_2", "P_2")
METHOD_LABELS <- c(
  "Before" = "Before\nmatching",
  "Full_1" = "Covariate\nonly",
  "Full_2" = "Single\npenalty",
  "P_2"    = "Dual\npenalty"
)
METHOD_COLORS <- c(
  "Before" = "#F8766D",
  "Full_1" = "#7CAE00",
  "Full_2" = "#00BFC4",
  "P_2"    = "#C77CFF"
)
BAND_COLORS <- c("d = 1" = "#D73027", "d = 2" = "#FC8D59", "d >= 3" = "#FEE08B")

BASE_SIZE <- 20 # slide-friendly

# =============================================================================
# LOAD CACHED TABLES
# =============================================================================
mahal_rds <- file.path(OUT_DIR, "eval_df_mahal.rds")
comp_rds  <- file.path(OUT_DIR, "eval_comp_agg.rds")
net_rds   <- file.path(OUT_DIR, "eval_df_net_plot.rds")

if (!file.exists(mahal_rds) || !file.exists(comp_rds) ||
    !file.exists(net_rds)) {
  stop("Cached tables not found. Run run_matching_eval.R first to produce ",
       "eval_df_mahal.rds, eval_comp_agg.rds and eval_df_net_plot.rds in:\n  ",
       OUT_DIR)
}

df_mahal <- readRDS(mahal_rds) %>% filter(dep_index == DEP)
comp_agg <- readRDS(comp_rds)  %>% filter(dep_index == DEP)
df_net   <- readRDS(net_rds)   %>% filter(dep_index == DEP)

cat(sprintf("[slide_eval] Loaded cached tables for dep=%d\n", DEP))

# =============================================================================
# PLOTS
# =============================================================================
p_mahal <- ggplot(df_mahal, aes(x = method, y = mahal_dist, fill = method)) +
  geom_boxplot(outlier.size = 0.9, outlier.alpha = 0.4, width = 0.55) +
  scale_fill_manual(values = METHOD_COLORS, guide = "none") +
  scale_x_discrete(labels = METHOD_LABELS, name = "Matching method") +
  scale_y_continuous(name = "Average distance") +
  labs(title = "Covariate dissimilarity") +
  theme_bw(base_size = BASE_SIZE) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x      = element_text(size = BASE_SIZE - 4))

plot_df <- comp_agg %>%
  filter(method == "Before" | scope == "within") %>%
  mutate(
    method_scope = dplyr::case_when(
      method == "Before" ~ "Before",
      TRUE               ~ paste0(as.character(method), "_within")
    ),
    method_scope = factor(method_scope, levels = c(
      "Before",
      "Full_1_within",
      "Full_2_within",
      "P_2_within"
    ))
  )

method_scope_labels <- c(
  "Before"        = "Before\nmatching",
  "Full_1_within" = "Covariate\nonly",
  "Full_2_within" = "Single\npenalty",
  "P_2_within"    = "Dual\npenalty"
)

p_net <- ggplot(df_net, aes(x = method, y = avg, fill = method)) +
  geom_boxplot(outlier.size = 0.9, outlier.alpha = 0.4, width = 0.55) +
  scale_fill_manual(values = METHOD_COLORS, guide = "none") +
  scale_x_discrete(labels = METHOD_LABELS, name = "Matching method") +
  scale_y_continuous(name = "Average distance") +
  labs(title = "Network distance") +
  theme_bw(base_size = BASE_SIZE) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x      = element_text(size = BASE_SIZE - 4))
ggsave(file.path(PLOTS, 'cov_mahal.pdf'), p_mahal, width=8, height=7)
ggsave(file.path(PLOTS, 'net_avg.pdf'), p_net, width=8, height=7)
all_mahal <- readRDS(mahal_rds)
all_net <- readRDS(net_rds)
p_mahal_all <- (p_mahal %+% all_mahal) + facet_wrap(~dep_label, ncol=3) + theme_bw(base_size=16) + theme(legend.position='none')
p_net_all <- (p_net %+% all_net) + facet_wrap(~dep_label, ncol=3) + theme_bw(base_size=16) + theme(legend.position='none')
ggsave(file.path(PLOTS, 'matching_eval_all_combined.pdf'), p_mahal_all / p_net_all, width=15, height=10)
stopifnot(identical(p_mahal$data, df_mahal), identical(p_net$data, df_net))

})