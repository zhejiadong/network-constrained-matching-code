#!/usr/bin/env Rscript
# =============================================================================
# plot_power.R
#
# Figure 1 (main): 3-panel line plot — Dual penalty (P_2) matching only
#   Panels: Independent | Moderate | Strong
#   x: beta_z, y: rejection rate, lines: 4 inference methods
#
# Figure 2 (appendix): 3x3 grid — all matching schemes x dep levels
#   Rows: Covariate only | Single penalty | Dual penalty
#   Cols: Independent | Moderate | Strong
#   Lines: 4 inference methods per panel
#
# Input:  power_outputs/ri_power_results.csv  (from run_ri_power.R)
# Output: power_outputs/power_figure.pdf
#         power_outputs/power_figure_appendix.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

DATA_DIR    <- "data"
OUT_DIR     <- file.path(DATA_DIR, "power_outputs")
IN_CSV      <- file.path(OUT_DIR, "ri_power_results.csv")
OUT_PDF     <- file.path(OUT_DIR, "power_figure.pdf")
OUT_PDF_APP <- file.path(OUT_DIR, "power_figure_appendix.pdf")

ALPHA <- 0.05

DEP_LEVELS  <- c("Independent", "Moderate", "Strong")
DEP_LABELS  <- c("1" = "Independent", "2" = "Moderate", "3" = "Strong")

SC_LEVELS   <- c("Covariate only", "Single penalty", "Dual penalty")
SC_LABELS   <- c("Full_1" = "Covariate only",
                 "Full_2" = "Single penalty",
                 "P_2"    = "Dual penalty")

METHODS_SHOW  <- c("LM_Crude", "LM_Subclass", "RI_Naive", "RI_Decay")
METHOD_LABELS <- c(
  "LM_Crude"    = "LM Crude",
  "LM_Subclass" = "LM Subclass",
  "RI_Naive"    = "RI Naive",
  "RI_Decay"    = "RI Sens"
)
METHOD_COLORS <- c(
  "LM_Crude"    = "#999999",
  "LM_Subclass" = "#E69F00",
  "RI_Naive"    = "#56B4E9",
  "RI_Decay"    = "#D55E00"
)
METHOD_LINES <- c(
  "LM_Crude"    = "dotted",
  "LM_Subclass" = "dashed",
  "RI_Naive"    = "dotdash",
  "RI_Decay"    = "solid"
)

# shared theme elements
theme_power <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      strip.background = element_rect(fill = "grey92", colour = "grey60"),
      strip.text       = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.x      = element_text(hjust = 0.5),
      legend.position  = "right"
    )
}

# shared scales
scale_method_colour <- function() {
  scale_colour_manual(
    values = METHOD_COLORS,
    breaks = METHODS_SHOW,
    labels = METHOD_LABELS[METHODS_SHOW],
    name   = "Method"
  )
}
scale_method_linetype <- function() {
  scale_linetype_manual(
    values = METHOD_LINES,
    breaks = METHODS_SHOW,
    labels = METHOD_LABELS[METHODS_SHOW],
    name   = "Method"
  )
}

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
cat("[plot_power] Reading", IN_CSV, "\n")
df <- read.csv(IN_CSV, stringsAsFactors = FALSE)

# LM_Crude is identical across subclasses — deduplicate
df_crude <- df |>
  filter(method == "LM_Crude") |>
  distinct(dep_index, bz, sim_id, method, p_val)

df_matched <- df |>
  filter(method != "LM_Crude", method %in% METHODS_SHOW)

# ---------------------------------------------------------------------------
# Helper: summarise to rejection rates
# ---------------------------------------------------------------------------
summarise_power <- function(data) {
  data |>
    group_by(dep_index, bz, subclass, method) |>
    summarise(
      reject_rate = mean(p_val < ALPHA, na.rm = TRUE),
      n           = n(),
      .groups     = "drop"
    ) |>
    mutate(
      dep_label = factor(DEP_LABELS[as.character(dep_index)],
                         levels = DEP_LEVELS),
      sc_label  = factor(SC_LABELS[subclass], levels = SC_LEVELS),
      method    = factor(method, levels = METHODS_SHOW)
    )
}

# ---------------------------------------------------------------------------
# Figure 1 (main): P_2 only
# ---------------------------------------------------------------------------
tbl_main <- bind_rows(
  df_crude  |> mutate(subclass = "P_2"),
  df_matched |> filter(subclass == "P_2")
) |> summarise_power()

cat("\n--- Power summary (P_2 / Dual penalty) ---\n")
print(as.data.frame(tbl_main |>
  select(dep_label, bz, method, reject_rate, n)),
  digits = 3, row.names = FALSE)

bz_breaks <- sort(unique(tbl_main$bz))

p_main <- ggplot(tbl_main,
    aes(x = bz, y = reject_rate,
        colour = method, linetype = method, group = method)) +
  geom_hline(yintercept = ALPHA, linetype = "dashed",
             colour = "black", linewidth = 0.9) +
  annotate("text", x = -Inf, y = ALPHA,
           label = "alpha * ' = 0.05'", parse = TRUE,
           hjust = 0, vjust = -0.4, size = 3.2, colour = "black") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_method_colour() +
  scale_method_linetype() +
  scale_x_continuous(breaks = bz_breaks, name = expression(beta[z])) +
  scale_y_continuous(name = "Rejection rate", limits = c(0, 1)) +
  facet_wrap(~ dep_label, ncol = 3) +
  theme_power() +
  labs(
    title = expression("Power by " * beta[z] *
                         "  (Dual-penalty matching)")
  )

ggsave(OUT_PDF, plot = p_main, width = 11, height = 4.5)
cat(sprintf("[plot_power] Main figure -> %s\n", OUT_PDF))

# ---------------------------------------------------------------------------
# Figure 2 (appendix): all 3 matching schemes, 3x3 grid
# ---------------------------------------------------------------------------
tbl_app <- bind_rows(
  df_crude   |> mutate(subclass = "Full_1") |>
    bind_rows(df_crude |> mutate(subclass = "Full_2")) |>
    bind_rows(df_crude |> mutate(subclass = "P_2")),
  df_matched
) |> summarise_power()

cat("\n--- Power summary (all matching schemes) ---\n")
print(as.data.frame(tbl_app |>
  select(dep_label, sc_label, bz, method, reject_rate, n)),
  digits = 3, row.names = FALSE)

p_app <- ggplot(tbl_app,
    aes(x = bz, y = reject_rate,
        colour = method, linetype = method, group = method)) +
  geom_hline(yintercept = ALPHA, linetype = "dashed",
             colour = "black", linewidth = 0.9) +
  annotate("text", x = -Inf, y = ALPHA,
           label = "alpha * ' = 0.05'", parse = TRUE,
           hjust = 0, vjust = -0.4, size = 3.2, colour = "black") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_method_colour() +
  scale_method_linetype() +
  scale_x_continuous(breaks = bz_breaks, name = expression(beta[z])) +
  scale_y_continuous(name = "Rejection rate", limits = c(0, 1)) +
  facet_grid(sc_label ~ dep_label) +
  theme_power() +
  labs(
    title = expression("Power by " * beta[z] *
                         " — all matching schemes")
  )

ggsave(OUT_PDF_APP, plot = p_app, width = 11, height = 9)
cat(sprintf("[plot_power] Appendix figure -> %s\n", OUT_PDF_APP))

cat("[plot_power] Done.\n")





