source("helpers.R")
input <- readRDS(file.path("data", paste0(s$mode, ".rds")))
dir.create("results", recursive = TRUE, showWarnings = FALSE)
pvalues <- grids <- diagnostics <- computation <- list()
power_eta <- c(0, 0.01, 0.03)
power_rho <- c(0, 0.05, 0.10)

for (i in seq_len(nrow(input$cases))) {
  case <- input$cases[i, ]
  data <- input$datasets[[i]]
  key <- data.frame(replication = case$replication,
                    dependence = dependence[case$dep], effect = case$effect)
  pvalues[[length(pvalues) + 1L]] <- cbind(key, design = "Unmatched",
    procedure = "Unmatched linear regression",
    p_value = summary(lm(Y ~ Z, data))$coefficients["Z", "Pr(>|t|)"])
  methods <- if (s$mode == "null") names(design_labels) else "dual"
  matches <- list()
  for (method in methods) {
    started <- proc.time()[["elapsed"]]
    m <- match_design(data, input$distance, method, s$solver)
    elapsed <- proc.time()[["elapsed"]] - started
    matches[[method]] <- m
    computation[[length(computation) + 1L]] <- cbind(key,
      design = design_labels[method], as.data.frame(m$solver_info),
      elapsed_seconds = elapsed)
    unadjusted <- netmatchRI::RI_unadjusted(m, "Y")$result$p_value
    regression <- summary(lm(Y ~ Z + factor(subclass), m$data))$coefficients["Z", "Pr(>|t|)"]
    values <- c("Linear regression with matched-set indicators" = regression,
                "Unadjusted randomization inference" = unadjusted)
    if (s$mode == "null") {
      set.seed(20260929L + case$dep * 100000L + case$replication * 10L + match(method, names(design_labels)))
      oracle <- oracle_test(m, input$covariance[[case$dep]], s$mc)
      values <- c(values, "Oracle randomization inference" = oracle$normal_p)
      computation[[length(computation)]]$oracle_normal_p <- oracle$normal_p
      computation[[length(computation)]]$oracle_variance <- oracle$variance
      if (method == "dual") {
        grid <- netmatchRI::sensitivity_grid(m, "Y",
          eta = seq(0.01, 0.05, by = 0.01), rho = seq(0.05, 0.40, by = 0.05))$grid
        grids[[length(grids) + 1L]] <- cbind(key[rep(1L, nrow(grid)), ], grid)
      }
    } else {
      adjusted <- netmatchRI::RI_adjusted(m, "Y",
        eta = power_eta[case$dep], rho = power_rho[case$dep])$result$p_value
      values <- c(values, "Adjusted randomization inference" = adjusted)
    }
    pvalues[[length(pvalues) + 1L]] <- cbind(key[rep(1L, length(values)), ],
      design = design_labels[method], procedure = names(values), p_value = unname(values))
  }
  if (s$mode == "null") {
    metrics <- matching_metrics(data, input$distance, matches)
    diagnostics[[i]] <- cbind(key[rep(1L, nrow(metrics)), ], metrics)
  }
  cat(s$mode, dependence[case$dep], "effect", case$effect,
      "replication", case$replication, "complete\n")
}
result <- list(settings = s, pvalues = do.call(rbind, pvalues),
               grid = do.call(rbind, grids), diagnostics = do.call(rbind, diagnostics),
               computation = do.call(rbind, computation))
saveRDS(result, file.path("results", paste0(s$mode, "_analysis.rds")))
