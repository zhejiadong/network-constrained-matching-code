
DATA_DIR <- "data"
load(file.path(DATA_DIR, "net.RData"))

make_W_row_standardized <- function(Adj) {
  Adj <- as.matrix(Adj)
  n <- nrow(Adj)
  
  # (recommended) remove self loops if any
  diag(Adj) <- 0
  
  # row sums = degrees (for unweighted adjacency)
  rs <- rowSums(Adj)
  
  W <- matrix(0, n, n)
  idx <- which(rs > 0)
  W[idx, ] <- Adj[idx, ] / rs[idx]   # each non-isolated row sums to 1
  
  W
}
W <- make_W_row_standardized(Adj)
generate_V_W_step2 <- function(W, kappa, alpha, addI = 1) {
  n <- nrow(W)
  I <- diag(n)
  B <- alpha * I + kappa * W
  Sigma <- B %*% t(B) + addI * I   # PSD
  Sigma
}
kappa <- c(0, 0.7,0.9)
alpha <- c(0,0.3,0.1)
V.list <- list()
scale_to_unit_diag <- function(V) {
  D <- sqrt(diag(V))
  V / (D %o% D)
}

for (i in 1:3) {
  V_raw <- generate_V_W_step2(
    W,
    kappa = kappa[i],
    alpha = alpha[i],
    addI = 0.001
  )
  
  V.list[[i]] <- scale_to_unit_diag(V_raw)
}
V.list[[1]] <- diag(popn)

# Re-save net.RData with V.list appended
net_file <- file.path(DATA_DIR, "net.RData")
load(net_file)
save(net_dist, net_proximity, Adj, popn, n.sim, n_blocks, block_sizes,
     kappa, alpha, V.list,
     file = net_file)