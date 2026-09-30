# Run from fhs_application/; use 02_analysis.R 2 or 02_analysis.R 4 for the supplement.
args <- commandArgs(trailingOnly = TRUE)
kappa <- if (length(args)) as.integer(args[1]) else 3L
stopifnot(kappa %in% c(2L, 3L, 4L))

required_packages <- c("netmatchRI", "ggplot2", "highs")
missing <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing))
  stop("Missing package(s): ", paste(missing, collapse = ", "), call. = FALSE)

input <- "data/fhs_analysis_data.rds"
if (!file.exists(input))
  stop("Run 01_prepare_data.R first; missing: ", input, call. = FALSE)
prepared <- readRDS(input)
if (!is.list(prepared) || !all(c("cohort", "distance") %in% names(prepared)))
  stop("Prepared FHS object has the wrong schema.", call. = FALSE)

cohort <- prepared$cohort
distance <- prepared$distance
treatment <- "obesity4"
outcome <- "lvm_height5"
covariates <- c("age4", "female", "smoking4", "alcohol_any4", "sbp4", "lvm_height4")
required <- c("id", treatment, outcome, covariates)
if (!is.data.frame(cohort) || nrow(cohort) != 2347L ||
    !all(required %in% names(cohort)) ||
    !all(stats::complete.cases(cohort[, required, drop = FALSE])))
  stop("Expected 2347 complete participant records with the required fields.", call. = FALSE)
if (sum(cohort[[treatment]] == 1) != 371L || sum(cohort[[treatment]] == 0) != 1976L)
  stop("Expected 371 participants with obesity and 1976 without.", call. = FALSE)
ids <- as.character(cohort$id)
if (!is.matrix(distance) || !identical(dim(distance), c(2347L, 2347L)) ||
    !identical(rownames(distance), ids) || !identical(colnames(distance), ids))
  stop("The network-distance matrix is not exactly aligned to the cohort.", call. = FALSE)

design_data <- cohort[, c("id", treatment, covariates), drop = FALSE]
rownames(design_data) <- ids
weight_type <- "ns"

cat("Running covariate-only full matching.\n")
covariate_match <- netmatchRI::netmatch(
  data = design_data,
  treat = treatment,
  covariates = covariates,
  network = distance,
  network_type = "distance",
  method = "covariate",
  kappa = kappa,
  min_controls = 0.01,
  max_controls = 100
)

cat("Running dual-penalty matching with kappa =", kappa, "\n")
dual_match <- netmatchRI::netmatch(
  data = design_data,
  treat = treatment,
  covariates = covariates,
  network = distance,
  network_type = "distance",
  method = "dual",
  kappa = kappa,
  solver = "highs",
  min_controls = 0.01,
  max_controls = 100,
  caliper = Inf,
  timelimit = 86400,
  mipgap = 0.01,
  threads = 1L
)
if (!inherits(covariate_match, "netmatch") || !inherits(dual_match, "netmatch"))
  stop("Matching did not return the expected netmatch objects.", call. = FALSE)

attach_outcome <- function(object) {
  index <- match(as.character(object$data$id), ids)
  if (anyNA(index)) stop("Could not align the outcome to a matched design.")
  object$data[[outcome]] <- cohort[[outcome]][index]
  object
}
covariate_match_y <- attach_outcome(covariate_match)
dual_match_y <- attach_outcome(dual_match)

cat("Running randomization inference and sensitivity analysis.\n")
covariate_naive <- netmatchRI::RI_unadjusted(
  covariate_match_y, outcome, kappa = kappa, weight_type = weight_type
)
dual_naive <- netmatchRI::RI_unadjusted(
  dual_match_y, outcome, kappa = kappa, weight_type = weight_type
)
dual_design <- netmatchRI::RI_design(
  dual_match_y, outcome, kappa = kappa, weight_type = weight_type
)
sensitivity <- netmatchRI::sensitivity_grid(
  dual_match_y,
  outcome,
  eta = seq(0.01, 0.05, by = 0.01),
  rho = seq(0, 1, by = 0.05),
  kappa = kappa,
  weight_type = weight_type
)
critical <- netmatchRI::critical_sensitivity(
  dual_match_y,
  outcome,
  rho = seq(0, 1, by = 0.01),
  alpha = 0.05,
  kappa = kappa,
  weight_type = weight_type
)

covariate_diagnostics <- netmatchRI::diagnose_match(covariate_match)
dual_diagnostics <- netmatchRI::diagnose_match(dual_match)
if (!is.data.frame(sensitivity$grid) || !is.data.frame(critical$curve))
  stop("Sensitivity analysis returned an unexpected object.", call. = FALSE)

extract_result <- function(x, design, inference) {
  out <- x$result
  if (!is.data.frame(out)) stop("RI result table is missing.")
  out$design <- design
  out$inference <- inference
  out
}
ri_results <- rbind(
  extract_result(covariate_naive, "Covariate-only", "Unadjusted RI"),
  extract_result(dual_naive, "Dual-penalty", "Unadjusted RI"),
  extract_result(dual_design, "Dual-penalty", "RI Design")
)

designs <- list("Covariate-only" = covariate_match, "Dual-penalty" = dual_match)
p_balance <- netmatchRI::plot_covariate_balance(designs)
p_similarity <- netmatchRI::plot_covariate_similarity(designs)
similarity <- p_similarity$data[, c("design", "average_mahalanobis")]

results_dir <- "results"
figures_dir <- "figures"
prefix <- paste0("fhs_e4_e5_k", kappa)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
results <- list(
  prepared_summary = prepared[c("cohort_flow", "graph_summary", "table1", "moran_report", "specification")],
  covariate_match = covariate_match,
  dual_match = dual_match,
  covariate_diagnostics = covariate_diagnostics,
  dual_diagnostics = dual_diagnostics,
  ri_results = ri_results,
  sensitivity = sensitivity,
  critical = critical,
  similarity = similarity,
  created = format(Sys.time(), tz = "UTC")
)
saveRDS(results, file.path(results_dir, paste0(prefix, "_results.rds")), compress = TRUE)
inference_summary <- data.frame(
  kappa = kappa, n = nrow(dual_match$data),
  n_sets = critical$detail$n_sets,
  covariate_unadjusted_p = covariate_naive$result$p_value,
  dual_unadjusted_p = dual_naive$result$p_value,
  design_p = dual_design$result$p_value,
  critical_eta_rho1 = critical$summary$eta_critical,
  relevant_set_pairs = critical$detail$n_relevant_pairs,
  critical_ratio = critical$summary$critical_ratio,
  critical_ratio_percent = 100 * critical$summary$critical_ratio,
  as.data.frame(dual_match$solver_info)
)
write.csv(inference_summary, file.path(results_dir, paste0(prefix, "_summary.csv")), row.names = FALSE)
print(inference_summary, row.names = FALSE)

p_sensitivity <- netmatchRI::plot_sensitivity(
  sensitivity, unadjusted = covariate_naive
)
ggplot2::ggsave(file.path(figures_dir, paste0(prefix, "_sensitivity_pvalues.pdf")),
                p_sensitivity, width = 7, height = 4)
ggplot2::ggsave(file.path(figures_dir, paste0(prefix, "_covariate_balance.pdf")),
                p_balance, width = 6.5, height = 4.5)
ggplot2::ggsave(file.path(figures_dir, paste0(prefix, "_covariate_similarity.pdf")),
                p_similarity, width = 6.5, height = 4.5)

outputs <- file.path(figures_dir, paste0(prefix, c(
  "_sensitivity_pvalues.pdf", "_covariate_balance.pdf", "_covariate_similarity.pdf"
)))
if (!all(file.exists(outputs)) || any(file.info(outputs)$size == 0))
  stop("One or more FHS figures were not created.", call. = FALSE)
cat("Saved FHS inference summary and three figures for kappa =", kappa, "\n")
