# Simulation settings: edit these values before running the three scripts.
s <- list(
  mode = "null",  # Change to "power" for the power analysis.
  reps = 500, # replication times
  n = 300, # sample size
  mc = 500, #Monte Carlo sampling times
  solver = "gurobi" # change to highs if needed
)

dependence <- c("Independent", "Moderate", "Strong")
design_labels <- c(covariate = "Covariate-only", single = "Single-penalty",
                   dual = "Dual-penalty")
covariates <- c("X1", "X2", "X3")

simulate_data <- function(V, effect, seed) {
  set.seed(seed)
  n <- nrow(V)
  X <- cbind(as.numeric(mvtnorm::rmvnorm(1, sigma = V)),
             as.numeric(mvtnorm::rmvnorm(1, sigma = V)),
             as.numeric(mvtnorm::rmvnorm(1, sigma = diag(n))))
  mu <- as.numeric(X %*% rep(0.1, 3))
  threshold <- uniroot(function(c) mean(1 - pnorm((c - mu) / sqrt(diag(V)))) - 0.23,
                       c(-60, 60))$root
  for (attempt in seq_len(100)) {
    Z <- as.integer(mvtnorm::rmvnorm(1, mean = mu, sigma = V) > threshold)
    if (sum(Z) >= 0.2 * n && sum(Z) <= 0.3 * n) break
  }
  Y <- as.numeric(mvtnorm::rmvnorm(1, mean = mu + effect * Z, sigma = V))
  data.frame(Y = Y, Z = Z, X1 = X[, 1], X2 = X[, 2], X3 = X[, 3])
}

match_design <- function(data, distance, method, solver) {
  run <- function(caliper = Inf, timelimit = 50, mipgap = 0.01) {
    netmatchRI::netmatch(data, "Z", covariates, distance,
      method = method, kappa = 2, solver = solver, network_type = "distance",
      min_controls = 0.01, max_controls = 100, caliper = caliper,
      timelimit = timelimit, mipgap = mipgap, threads = 1)
  }
  if (method != "dual") return(run())
  # The second pass relaxes the covariate caliper, not the network constraint.
  tryCatch(run(8, 50, 0.01), error = function(e) {
    message("Retrying dual matching: ", conditionMessage(e))
    run(12, 150, 0.02)
  })
}

oracle_test <- function(match, V, mc = 500) {
  df <- match$data
  ids <- as.integer(rownames(df))
  sets <- split(seq_len(nrow(df)), df$subclass, drop = TRUE)
  score <- numeric(nrow(df))
  constant <- expectation <- 0
  for (i in sets) {
    nt <- sum(df$Z[i])
    w <- 1 / (length(i) + 1)
    score[i] <- w * rank(df$Y[i])
    constant <- constant + w * nt * (nt + 1) / 2
    expectation <- expectation + w * nt * (length(i) - nt) / 2
  }
  observed <- sum(df$Z * score) - constant
  G <- crossprod(chol(V[ids, ids, drop = FALSE]),
                 matrix(rnorm(length(ids) * mc), length(ids), mc))
  values <- vapply(seq_len(mc), function(b) {
    z <- integer(length(ids))
    for (i in sets) {
      treated <- i[order(G[i, b], decreasing = TRUE)[seq_len(sum(df$Z[i]))]]
      z[treated] <- 1L
    }
    sum(z * score) - constant
  }, numeric(1))
  list(p_value = (1 + sum(abs(values - expectation) >= abs(observed - expectation))) / (mc + 1),
       normal_p = 2 * pnorm(-abs((observed - expectation) / sqrt(var(values)))),
       variance = var(values), values = values)
}

matching_metrics <- function(data, distance, matches) {
  X <- as.matrix(data[, covariates])
  S <- cov(X)
  Sinv <- if (abs(det(S)) > 1e-10) solve(S) else MASS::ginv(S)
  Q <- X %*% Sinv %*% t(X)
  mahal <- sqrt(pmax(outer(diag(Q), diag(Q), "+") - 2 * Q, 0))
  z <- data$Z
  pooled_sd <- vapply(covariates, function(v)
    sqrt((var(data[z == 1, v]) + var(data[z == 0, v])) / 2), numeric(1))
  rows <- list()
  for (method in c("Before", names(matches))) {
    s <- if (method == "Before") rep(1L, nrow(data)) else matches[[method]]$subclass
    same <- outer(as.character(s), as.character(s), "==")
    tc <- outer(z, z, function(a, b) a == 1 & b == 0) & same
    nd <- distance[upper.tri(distance) & same]
    sets <- split(seq_len(nrow(data)), s, drop = TRUE)
    smd <- vapply(seq_along(covariates), function(j) {
      differences <- vapply(sets, function(i)
        mean(X[i[z[i] == 1], j]) - mean(X[i[z[i] == 0], j]), numeric(1))
      abs(weighted.mean(differences, lengths(sets)) / pooled_sd[j])
    }, numeric(1))
    rows[[method]] <- data.frame(design = if (method == "Before") method else design_labels[method],
      metric = c("Covariate dissimilarity", "Network distance", paste0("SMD: ", covariates)),
      value = c(mean(mahal[tc]), mean(nd[is.finite(nd)]), smd), row.names = NULL)
  }
  do.call(rbind, rows)
}
