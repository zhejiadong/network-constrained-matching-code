# Re-run the manuscript's null inference on its frozen original matched data.
# Original helper implementations are in source/, unchanged. No matching occurs.
args <- commandArgs(trailingOnly=TRUE)
root <- '.'
if (!length(args)) stop('Usage: Rscript 06_run_null_inference.R TASK [smoke]')
task <- as.integer(args[1]); smoke <- length(args)>1 && args[2]=='smoke'
stopifnot(task %in% 1:150)
dep <- (task-1L)%/%50L+1L
first <- ((task-1L)%%50L)*10L+1L
reps <- if(smoke) first else first:(first+9L)
input_file <- file.path(root,'inputs',sprintf('null_dep%d.rds',dep))
input <- readRDS(input_file)
source(file.path(root,'source/statistics.R'))
source(file.path(root,'source/covariance.R'))
source(file.path(root,'source/utility_oracle.R'))
# Cache the original deterministic 10,000-quantile calculation by set counts.
# This changes repeated work only, not the historical numerical formula.
uncached_bound <- get_design_cov_bound
cache <- new.env(parent=emptyenv())
get_design_cov_bound <- function(nt_s,nc_s,nt_k,nc_k) {
 key<-paste(nt_s,nc_s,nt_k,nc_k,sep='_')
 if(!exists(key,cache,inherits=FALSE)) assign(key,uncached_bound(nt_s,nc_s,nt_k,nc_k),cache)
 get(key,cache,inherits=FALSE)
}
source_files <- c(file.path(root,'06_run_null_inference.R'),
                  file.path(root,'source',c('statistics.R','covariance.R','utility_oracle.R')))
identity <- list(input_md5=unname(tools::md5sum(input_file)),source_md5=tools::md5sum(source_files),
                 task=task,dep=dep,reps=reps,mc=500L,seed_base=20260929L,scope='original-matched-null')
out_dir<-file.path(root,if(smoke)'smoke' else 'tasks');dir.create(out_dir,showWarnings=FALSE)
out_file<-file.path(out_dir,sprintf('task_%03d.rds',task))
if(file.exists(out_file)) {
 old<-readRDS(out_file);stopifnot(identical(old$identity,identity))
 cat('PASS existing checkpoint',out_file,'\n');quit(save='no')
}
set_cols<-c(Full_1='subclass_full1',Full_2='subclass_full2',P_2='subclass_p_2')
tables<-grids<-list(); started<-Sys.time()
for(i in reps) {
 df<-input$datasets[[i]]
 stopifnot(nrow(df)==300L,identical(rownames(df),as.character(1:300)))
 crude<-summary(lm(Y~Z,data=df))$coefficients['Z','Pr(>|t|)']
 for(k in seq_along(set_cols)) {
  design<-names(set_cols)[k]; x<-df;x$subclass<-df[[set_cols[k]]]
  stopifnot(!anyNA(x$subclass),!anyDuplicated(x$Y))
  obs<-get_observed_stats(x,'ns','subclass','Z','Y')
  st<-obs$set_stats;stopifnot(all(st$nt>0),all(st$nc>0))
  diag_var<-sum(st$weight^2*st$nt*st$nc*(st$nt+st$nc+1)/12)
  naive<-2*(1-pnorm(abs((obs$T_obs-obs$E_T)/sqrt(diag_var))))
  lm_sub<-summary(lm(Y~Z+factor(subclass),data=x))$coefficients['Z','Pr(>|t|)']
  # Original implementation, B=500. A new explicit seed makes this rerun resumable.
  seed<-20260929L+dep*100000L+i*10L+k;set.seed(seed)
  mc<-mc_var_T_oracle(x,input$covariance,mc=500L)
  oracle_z<-(obs$T_obs-obs$E_T)/sqrt(mc$var_T)
  p_norm<-2*(1-pnorm(abs(oracle_z)))
  p_emp<-(1+sum(abs(mc$T_vals-obs$E_T)>=abs(obs$T_obs-obs$E_T)))/501
  tables[[length(tables)+1L]]<-data.frame(dep_index=dep,sim_id=i,design=design,
   LM_Crude=crude,LM_Subclass=lm_sub,RI_Naive=naive,Oracle_Normal=p_norm,Oracle_Empirical=p_emp,
   T_obs=obs$T_obs,E_T=obs$E_T,oracle_variance=mc$var_T,oracle_seed=seed)
  if(design=='P_2') {
   M<-precompute_design_cov_matrix(x,'subclass','Z')
   D<-get_set_distance_matrix(input$network_distance,x$subclass,'dmin')
   grid<-expand.grid(eta=seq(.01,.05,.01),rho=seq(.05,.40,.05))
   grid$Variance<-vapply(seq_len(nrow(grid)),function(j) {
    W<-get_set_decay_matrix(D,grid$eta[j],grid$rho[j],2,diag_val=1)
    as.numeric(t(st$weight)%*%get_U_cov_from_precomputed(M,W)%*%st$weight)
   },numeric(1))
   grid$p_value<-2*(1-pnorm(abs((obs$T_obs-obs$E_T)/sqrt(grid$Variance))))
   grid$dep_index<-dep;grid$sim_id<-i;grid$design<-design
   grids[[length(grids)+1L]]<-grid
  }
 }
 cat(sprintf('task=%d dep=%d replication=%d complete\n',task,dep,i));flush.console()
}
tab<-do.call(rbind,tables);grid<-do.call(rbind,grids)
stopifnot(nrow(tab)==length(reps)*3L,nrow(grid)==length(reps)*40L,
 all(vapply(tab[vapply(tab,is.numeric,logical(1))],function(x)all(is.finite(x)),logical(1))),all(is.finite(grid$p_value)))
result<-list(identity=identity,table=tab,grid=grid,elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')))
tmp<-paste0(out_file,'.tmp');saveRDS(result,tmp);stopifnot(file.rename(tmp,out_file))
cat('PASS:',out_file,'elapsed_seconds=',result$elapsed_seconds,'\n')
