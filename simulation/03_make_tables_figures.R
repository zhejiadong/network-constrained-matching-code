library(ggplot2)
source("helpers.R")
result <- readRDS(file.path("results", paste0(s$mode, "_analysis.rds")))
results_dir <- "results"
figures_dir <- "figures"
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
rows <- result$pvalues
rows$reject <- rows$p_value < 0.05
rates <- aggregate(reject ~ dependence + effect + design + procedure, rows, mean)
n <- aggregate(reject ~ dependence + effect + design + procedure, rows, length)
rates$n <- n$reject
names(rates)[names(rates) == "reject"] <- "rate"
rates$dependence <- factor(rates$dependence, levels = dependence)
write.csv(rates, file.path(results_dir, paste0(s$mode, "_rejection_rates.csv")), row.names = FALSE)
write.csv(result$computation, file.path(results_dir, paste0(s$mode, "_computation.csv")), row.names = FALSE)
theme_set(theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank()))

if (s$mode == "power") {
  p <- ggplot(rates, aes(effect, rate, colour = procedure, linetype = procedure)) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    geom_line() + geom_point() + facet_wrap(~ dependence, nrow = 1) +
    scale_x_continuous(breaks = seq(0, 0.4, 0.1)) +
    labs(x = expression(beta[z]), y = "Rejection rate", colour = NULL, linetype = NULL) +
    theme(legend.position = "bottom")
  ggsave(file.path(figures_dir, "power_figure.pdf"), p, width = 11, height = 4.5)
} else {
  table <- reshape(rates[, c("design", "procedure", "dependence", "rate")],
    idvar = c("design", "procedure"), timevar = "dependence", direction = "wide")
  names(table) <- sub("^rate[.]", "", names(table))
  table <- table[, c("design", "procedure", dependence)]
  write.csv(table, file.path(results_dir, "null_table.csv"), row.names = FALSE)
  tex <- c("\\begin{tabular}{llrrr}", "\\hline",
           "Matching design & Inference procedure & Independent & Moderate & Strong \\\\", "\\hline")
  for (i in seq_len(nrow(table)))
    tex <- c(tex, paste0(paste(c(table$design[i], table$procedure[i],
      sprintf("%.3f", as.numeric(table[i, dependence]))), collapse = " & "), " \\\\"))
  writeLines(c(tex, "\\hline", "\\end{tabular}"), file.path(results_dir, "null_table.tex"))

  grid <- result$grid
  grid$reject <- grid$p_value < 0.05
  grid <- aggregate(reject ~ dependence + eta + rho, grid, mean)
  grid$dependence <- factor(grid$dependence, levels = dependence)
  baseline <- rates[rates$design == "Dual-penalty" &
                      rates$procedure == "Unadjusted randomization inference", ]
  p <- ggplot(grid, aes(rho, reject, colour = factor(eta), group = eta)) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    geom_hline(data = baseline, aes(yintercept = rate), colour = "red3", linetype = "dotted") +
    geom_line() + geom_point() + facet_wrap(~ dependence, nrow = 1) +
    labs(x = expression(rho), y = "Rejection rate", colour = expression(eta))
  ggsave(file.path(figures_dir, "ri_grid_linePlot.pdf"), p, width = 11, height = 4)
  write.csv(grid, file.path(results_dir, "sensitivity_rejection_rates.csv"), row.names = FALSE)

  d <- result$diagnostics
  d$dependence <- factor(d$dependence, levels = dependence)
  d$design <- factor(d$design, levels = c("Before", unname(design_labels)))
  boxplot_metric <- function(data, metric) {
    ggplot(data[data$metric == metric, ], aes(design, value, fill = design)) +
      geom_boxplot(width = 0.6) + labs(x = NULL, y = metric) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
  }
  strong <- d[d$dependence == "Strong", ]
  ggsave(file.path(figures_dir, "cov_mahal.pdf"),
    boxplot_metric(strong, "Covariate dissimilarity"), width = 5, height = 4)
  ggsave(file.path(figures_dir, "net_avg.pdf"),
    boxplot_metric(strong, "Network distance"), width = 5, height = 4)
  plotted <- d[d$metric %in% c("Covariate dissimilarity", "Network distance"), ]
  p <- ggplot(plotted, aes(design, value, fill = design)) + geom_boxplot(width = 0.6) +
    facet_grid(metric ~ dependence, scales = "free_y") + labs(x = NULL, y = NULL) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(figures_dir, "matching_eval_all_combined.pdf"), p, width = 11, height = 7)
  write.csv(d, file.path(results_dir, "matching_diagnostics.csv"), row.names = FALSE)
}
cat("Saved", s$mode, "tables and figures.\n")
