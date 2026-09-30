# Run from simulation/ after matching, or use accepted saved null inputs directly.
source_root <- 'data'; out <- 'inputs'
dir.create(out, showWarnings=FALSE)
net <- new.env(); load(file.path(source_root,'net.RData'),net)
for (d in 1:3) {
 e <- new.env();load(file.path(source_root,sprintf('matched_multiprox_bz_0_v%d.RData',d)),e)
 stopifnot(length(e$df.final.list)==500L)
 for (x in e$df.final.list) stopifnot(nrow(x)==300L, identical(rownames(x),as.character(1:300)), !anyNA(x), all(x$Z %in% 0:1))
 x <- list(dep_index=d, datasets=e$df.final.list, covariance=net$V.list[[d]], network_distance=net$net_dist)
 dest<-file.path(out,sprintf('null_dep%d.rds',d));saveRDS(x,dest)
 stopifnot(identical(x,readRDS(dest)));cat('PASS:',dest,'\n')
}
