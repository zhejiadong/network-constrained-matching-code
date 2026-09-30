library(tidyr)
library(data.table)
library(dplyr)
library(energy)
library(sna)
library(mvtnorm)
library(RItools)
library(igraph)
library(ape)
library(Matrix)

DATA_DIR <- "data"

# Load network and covariance objects saved by net.R and variance.R
load(file.path(DATA_DIR, "net.RData"))


alpha_vec <- c(1,1,0)
#beta.x =  matrix(c(0.1, 0.1, 0.1,0.1, 0.1), nrow = 5, ncol = 1) # X -> Y
beta.x =  matrix(c(0.1, 0.1, 0.1), nrow = 3, ncol = 1) # X -> Y
beta.z = 0 # Z -> Y
gamma.x <- matrix(c(0.1,0.1,0.1), nrow = 3, ncol = 1) # X -> Z
#gamma.x <- matrix(c(0.1,0.1,0.1,0.1,0.1), nrow = 5, ncol = 1) # X -> Z
sim.process <- function(V,beta.x, gamma.x, beta.z,alpha_vec,
                        seeds = NULL,
                        target.treat = FALSE,K_min, K_max, max_try = 1000){
  
  popn <- nrow(V)
  
  if (!is.null(seeds)) set.seed(seeds)
  X <- matrix(0, nrow = popn, ncol = length(alpha_vec))  
  for(j in 1: length(alpha_vec)){
    alpha <- alpha_vec[j]
    Sigma_X <- alpha * V + (1 - alpha) * diag(popn)
    X[,j] <- rmvnorm(1, mean = rep(0, popn), sigma = Sigma_X)
  }
  attempt <- 1
  
  mu <- as.numeric(X %*% gamma.x)
  
  
  if(!isFALSE(target.treat)){
    c_thr <- solve_c_global(mu, V, target.treat)
  } 
  repeat {
    Z_lin <- as.numeric(rmvnorm(1, mean = mu, sigma = V))
    
    if(!isFALSE(target.treat)){
      Z <- as.integer(Z_lin > c_thr)
    } else {
      Z <- as.integer(Z_lin > 0)
      c_thr <- 0
    }
    n1 <- sum(Z)
    if ((n1 >= K_min && n1 <= K_max) || attempt >= max_try) break
    attempt <- attempt + 1
  }
  if (n1 < K_min || n1 > K_max)
    warning(sprintf("Reached %d attempts but only %d treated units (min=%d)", 
                    attempt, n1, K_min))
  
  
  Y <- rmvnorm(1, mean = X %*% beta.x + Z * beta.z, sigma = V) %>% as.numeric()
  res <- cbind(Y,Z,X) %>% data.frame()
  names(res) <- c("Y","Z",paste0("X",1:length(alpha_vec)))
  return(list(
    "df" = res
  ))
}

solve_c_global <- function(mu, V, p) {
  sd_i <- sqrt(diag(V))
  g <- function(c) mean(1 - pnorm((c - mu) / sd_i)) - p
  uniroot(g, c(-60, 60))$root   
}

beta.z.Vec <- c(0,0.1,0.2,0.3,0.4)
# This packaged entrypoint selects the accepted null-only scope.
for(b in 1L){
  beta.z <- beta.z.Vec[b]
  for(dep_indx in c(1:3)){
    V <- V.list[[dep_indx]]
    
    df.list <- list()
    lm.est <- vector()
    lm.pvalue <- vector()
    MoranZ <- vector()
    for(i in 1:n.sim){
      print(paste0(beta.z,":",dep_indx,":",i))
      seeds <- 123 + i + beta.z *1e3 + kappa[dep_indx] * 1e5
      per.sim  <-  sim.process(V, beta.x, gamma.x, 
                               alpha_vec = alpha_vec,
                               beta.z,seeds = seeds,
                               target.treat = 0.23, K_min = 0.2*popn, K_max = 0.3*popn,max_try = 100)
      df.list[[i]] <- per.sim$df
      #lm.model <- lm(Y~ Z+ X1 + X2 + X3 + X4 + X5, data = df.list[[i]])
      lm.model <- lm(Y~ Z+ X1 + X2 + X3, data = df.list[[i]])
      lm.est[i]<- lm.model$coefficients[["Z"]]
      s.model <- summary(lm.model)
      lm.pvalue[i] <- s.model$coefficients[2,4]
      MoranZ[i] <- Moran.I(df.list[[i]]$Z, Adj)$p.value
    }

    kappa.true = kappa[dep_indx]
    alpha.true = alpha[dep_indx]
    save( df.list,
          lm.est, n.sim, popn, 
          lm.pvalue,
          moranZ = MoranZ,
          net_dist, net_proximity,
          dep_indx = dep_indx,
          kappa.true = kappa.true,
          alpha.true = alpha.true,
          n_blocks,
          block_sizes,
          beta.z,
          beta.x,
          gamma.x,
          file = file.path(DATA_DIR, paste0("bz", beta.z, "_v", dep_indx, ".RData")))
    
  }
}


