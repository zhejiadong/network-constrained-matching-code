# Direct analysis of the locally prepared FHS application object.
# Run 01_prepare_data.R first. All generated objects remain in ignored local folders.

script_path <- function() {
  frame_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(frame_file)) return(normalizePath(frame_file, winslash = "/"))
  arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(arg)) return(normalizePath(sub("^--file=", "", arg[1]), winslash = "/"))
  normalizePath("02_analysis.R", winslash = "/", mustWork = TRUE)
}

required_packages <- c("netmatchRI", "ggplot2", "MASS", "highs")
missing <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing))
  stop("Missing package(s): ", paste(missing, collapse = ", "), call. = FALSE)

base <- dirname(script_path())
input <- file.path(base, "data", "fhs_analysis_data.rds")
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
  stop("The prepared cohort does not satisfy the reviewed schema/N contract.", call. = FALSE)
if (sum(cohort[[treatment]] == 1) != 371L || sum(cohort[[treatment]] == 0) != 1976L)
  stop("Treatment counts differ from the reviewed 371/1,976 contract.", call. = FALSE)
ids <- as.character(cohort$id)
if (!is.matrix(distance) || !identical(dim(distance), c(2347L, 2347L)) ||
    !identical(rownames(distance), ids) || !identical(colnames(distance), ids))
  stop("The network-distance matrix is not exactly aligned to the cohort.", call. = FALSE)

design_data <- cohort[, c("id", treatment, covariates), drop = FALSE]
rownames(design_data) <- ids
kappa <- 3
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

cat("Running dual-penalty matching (kappa = 3).\n")
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
  timelimit = 7200,
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
  eta = c(0, 0.02, 0.05, 0.10, 0.15, 0.20),
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
  extract_result(covariate_naive, "Covariate-only", "RI Naive"),
  extract_result(dual_naive, "Dual-penalty", "RI Naive"),
  extract_result(dual_design, "Dual-penalty", "RI Design")
)

inverse_covariance <- function(X) {
  S <- stats::cov(X)
  if (abs(det(S)) > 1e-10) solve(S) else MASS::ginv(S)
}
average_mahalanobis <- function(object = NULL) {
  X <- as.matrix(cohort[, covariates, drop = FALSE])
  Sinv <- inverse_covariance(X)
  z <- cohort[[treatment]]
  if (is.null(object)) {
    pairs <- expand.grid(treated = which(z == 1), control = which(z == 0))
  } else {
    index <- match(object$data$id, ids)
    subclass <- object$data$subclass
    sets <- split(seq_along(index), subclass)
    pairs <- do.call(rbind, lapply(sets, function(position) {
      units <- index[position]
      tr <- units[z[units] == 1]
      co <- units[z[units] == 0]
      if (!length(tr) || !length(co)) return(NULL)
      expand.grid(treated = tr, control = co)
    }))
  }
  values <- vapply(seq_len(nrow(pairs)), function(i) {
    delta <- X[pairs$treated[i], ] - X[pairs$control[i], ]
    sqrt(as.numeric(t(delta) %*% Sinv %*% delta))
  }, numeric(1))
  mean(values)
}
similarity <- data.frame(
  design = c("Before", "Covariate-only", "Dual-penalty"),
  average_mahalanobis = c(
    average_mahalanobis(),
    average_mahalanobis(covariate_match),
    average_mahalanobis(dual_match)
  ),
  stringsAsFactors = FALSE
)

results_dir <- file.path(base, "results")
figures_dir <- file.path(base, "figures")
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
saveRDS(results, file.path(results_dir, "fhs_analysis_results.rds"), compress = TRUE)

base_theme <- ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

grid <- sensitivity$grid
grid$eta_label <- factor(sprintf("%.2f", grid$eta),
  levels = sprintf("%.2f", sort(unique(grid$eta))))
naive_p <- ri_results$p_value[
  ri_results$design == "Covariate-only" & ri_results$inference == "RI Naive"
]
p_sensitivity <- ggplot2::ggplot(grid,
  ggplot2::aes(rho, p_value, colour = eta_label, group = eta_label)) +
  ggplot2::geom_line(linewidth = 0.65) + ggplot2::geom_point(size = 1.6) +
  ggplot2::geom_hline(yintercept = 0.05, colour = "red", linetype = "dotted") +
  ggplot2::geom_hline(yintercept = naive_p, linetype = "dashed") +
  ggplot2::labs(x = expression(rho), y = "p-value", colour = expression(eta)) + base_theme
ggplot2::ggsave(file.path(figures_dir, "fhs_sensitivity_pvalues.pdf"),
                p_sensitivity, width = 7, height = 4)

curve <- critical$curve
curve$eta_display <- pmin(1, pmax(0, curve$eta_critical))
p_critical <- ggplot2::ggplot(curve, ggplot2::aes(rho, eta_display)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = eta_display),
                       fill = "#D9EAF7", alpha = 0.65) +
  ggplot2::geom_line(colour = "#2166AC", linewidth = 0.7) +
  ggplot2::geom_point(data = data.frame(rho = 1, eta_display = 1),
                      shape = 18, size = 2.6) +
  ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(x = expression(rho), y = expression(eta^"*")) + base_theme
ggplot2::ggsave(file.path(figures_dir, "fhs_critical_sensitivity.pdf"),
                p_critical, width = 7, height = 4)

cov_balance <- covariate_diagnostics$covariate_balance
dual_balance <- dual_diagnostics$covariate_balance
balance <- rbind(
  data.frame(covariate = cov_balance$covariate, design = "Before",
             abs_smd = cov_balance$before_abs_smd),
  data.frame(covariate = cov_balance$covariate, design = "Covariate-only",
             abs_smd = cov_balance$after_abs_smd),
  data.frame(covariate = dual_balance$covariate, design = "Dual-penalty",
             abs_smd = dual_balance$after_abs_smd)
)
balance$design <- factor(balance$design,
  levels = c("Before", "Covariate-only", "Dual-penalty"))
p_balance <- ggplot2::ggplot(balance,
  ggplot2::aes(abs_smd, factor(covariate), colour = design, shape = design)) +
  ggplot2::geom_vline(xintercept = 0.1, linetype = "dashed") +
  ggplot2::geom_point(size = 2.6) +
  ggplot2::labs(x = "Absolute standardized mean difference", y = NULL,
                colour = NULL, shape = NULL) + base_theme +
  ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(figures_dir, "fhs_covariate_balance.pdf"),
                p_balance, width = 6.5, height = 4.5)

similarity$design <- factor(similarity$design,
  levels = c("Before", "Covariate-only", "Dual-penalty"))
p_similarity <- ggplot2::ggplot(similarity,
  ggplot2::aes(design, average_mahalanobis, fill = design)) +
  ggplot2::geom_col(width = 0.62, colour = "#333333") +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", average_mahalanobis)),
                     vjust = -0.45) +
  ggplot2::labs(x = NULL, y = "Average treated-control Mahalanobis distance") +
  base_theme + ggplot2::theme(legend.position = "none")
ggplot2::ggsave(file.path(figures_dir, "fhs_covariate_similarity.pdf"),
                p_similarity, width = 6.5, height = 4.5)

outputs <- file.path(figures_dir, c(
  "fhs_sensitivity_pvalues.pdf", "fhs_critical_sensitivity.pdf",
  "fhs_covariate_balance.pdf", "fhs_covariate_similarity.pdf"
))
if (!all(file.exists(outputs)) || any(file.info(outputs)$size == 0))
  stop("One or more FHS figures were not created.", call. = FALSE)
cat("PASS: FHS matching, inference, sensitivity, and four figures completed.\n")
