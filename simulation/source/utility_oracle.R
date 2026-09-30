draw_Z_oracle<- function(df, V, s_col = "subclass", z_obs_col = "Z", mu = NULL, g = NULL) {
  n <- nrow(df)
  if (is.null(mu)) mu <- rep(0, n)
  if (is.null(g)) {
    g <- as.numeric(rmvnorm(1, mean = mu, sigma = V))
  }
  sets <- sort(unique(df[[s_col]]))
  Z_new <- rep(0L, n)
  for (s in sets) {
    ids <- which(df[[s_col]] == s)
    nt <- sum(df[[z_obs_col]][ids] == 1)
    ns <- length(ids)
    if (nt <= 0 || nt >= ns) {
      if (nt >= ns) Z_new[ids] <- 1L
      next
    }
    g_sub <- g[ids]
    ord <- order(g_sub, decreasing = TRUE)
    tr <- ids[ord[seq_len(nt)]]
    Z_new[tr] <- 1L
  }
  Z_new
}

mc_var_T_oracle <- function(df, V, mc = 2000, subclass_col = "subclass",
                        treat_col = "Z", outcome_col = "Y",
                        weight_type = "ns", mu = NULL) {
  base_stats <- get_observed_stats(df, weight_type, subclass_col, treat_col, outcome_col)
  w <- base_stats$set_stats$weight
  n <- nrow(df)
  if (is.null(mu)) mu <- rep(0, n)

  # --- Suggestion 1: precompute per-unit rank scores (Y fixed under null) ---
  # T = sum_s w_s * U_s, and U_s = W_s - n_t*(n_t+1)/2
  # where W_s = sum of within-set Y ranks of treated units.
  # So T = sum_i Z_i * score_i - constant, computed in O(n) per MC draw.
  sets      <- sort(unique(df[[subclass_col]]))
  score_vec <- numeric(n)   # w_s * rank_{s,i} for each unit i
  constant  <- 0            # sum_s w_s * n_t_s * (n_t_s + 1) / 2
  for (k in seq_along(sets)) {
    ids  <- which(df[[subclass_col]] == sets[k])
    Y_s  <- df[[outcome_col]][ids]
    n_t  <- sum(df[[treat_col]][ids])
    w_s  <- w[k]
    score_vec[ids] <- w_s * rank(Y_s, ties.method = "average")
    constant <- constant + w_s * n_t * (n_t + 1) / 2
  }

  # --- Suggestion 2: generate all MC g-vectors at once (one BLAS call) ---
  L <- tryCatch(chol(V), error = function(e) NULL)
  use_chol <- !is.null(L)
  Z_raw <- matrix(rnorm(n * mc), n, mc)
  if (use_chol) {
    G_mat <- crossprod(L, Z_raw)          # t(L) %*% Z_raw, triangular BLAS
  } else {
    ev    <- eigen(V, symmetric = TRUE)
    sqrtV <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0))) %*% t(ev$vectors)
    G_mat <- sqrtV %*% Z_raw
  }
  G_mat <- G_mat + mu                     # add mean (no-op when mu = 0)

  # --- MC loop: only Z assignment + O(n) dot product ---
  T_vals <- numeric(mc)
  for (b in seq_len(mc)) {
    Z_new     <- draw_Z_oracle(df, V, subclass_col, treat_col, mu, g = G_mat[, b])
    T_vals[b] <- sum(Z_new * score_vec) - constant
  }

  list(T_vals = T_vals, var_T = var(T_vals), weights = w)
}
