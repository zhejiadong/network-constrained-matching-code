source("helpers.R")
dir.create("data", recursive = TRUE, showWarnings = FALSE)

set.seed(2026)
probabilities <- matrix(0.003, 4, 4)
diag(probabilities) <- 0.19
graph <- igraph::sample_sbm(s$n, pref.matrix = probabilities,
                          block.sizes = rep(s$n / 4, 4), directed = FALSE)
adjacency <- as.matrix(igraph::as_adjacency_matrix(graph))
distance <- igraph::distances(graph)
W <- adjacency / pmax(rowSums(adjacency), 1)
alpha1 <- c(0, 0.7, 0.9)
alpha2 <- c(1, 0.3, 0.1)
V <- lapply(seq_along(alpha1), function(j) {
  B <- alpha1[j] * W + alpha2[j] * diag(s$n)
  # Keep the original simulation's ridge and unit marginal variances.
  raw <- tcrossprod(B) + 0.001 * diag(s$n)
  raw / outer(sqrt(diag(raw)), sqrt(diag(raw)))
})
V[[1]] <- diag(s$n)
effects <- if (s$mode == "power") seq(0, 0.4, by = 0.1) else 0
cases <- expand.grid(replication = seq_len(s$reps), dep = 1:3, effect = effects)
datasets <- lapply(seq_len(nrow(cases)), function(i) {
  x <- cases[i, ]
  simulate_data(V[[x$dep]], x$effect,
                123 + x$replication + x$effect * 1000 + alpha1[x$dep] * 100000)
})
saveRDS(list(settings = s, cases = cases, datasets = datasets,
             distance = distance, covariance = V),
        file.path("data", paste0(s$mode, ".rds")))
cat("Generated", nrow(cases), "datasets on a fixed", s$n, "node network.\n")
