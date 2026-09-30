library(igraph)

DATA_DIR <- "data"
dir.create(DATA_DIR, showWarnings = FALSE)

popn <- 300
n.sim = 500
n_blocks <- 4
block_sizes <- rep(popn / n_blocks, n_blocks)

p_in  <- 0.19
p_out <- 0.003
pm <- matrix(p_out, nrow = n_blocks, ncol = n_blocks)
diag(pm) <- p_in

set.seed(2026)
g <- sample_sbm(popn, pref.matrix = pm, block.sizes = block_sizes, directed = FALSE)

pdf(file.path(DATA_DIR, "network.pdf"))
plot(g, vertex.label = NA, vertex.size = 3,
     cex.main = 1.5 ,cex.sub = 1.5
)
dev.off()
net_dist <- igraph::distances(g)
net_proximity <- 1/net_dist

igraph::degree(g) %>% mean()
Adj <- as.matrix(as_adjacency_matrix(g))

1/(eigen(Adj)$values %>% range())
mean(net_dist[ net_dist>0])
mean(net_dist[is.finite(net_dist) & net_dist>0])
table(net_dist[is.finite(net_dist) & net_dist>0])/sum(table(net_dist[is.finite(net_dist) & net_dist>0]))

save(net_dist, net_proximity, Adj, popn, n.sim, n_blocks, block_sizes,
     file = file.path(DATA_DIR, "net.RData"))
