# Summarize accepted null task outputs; no simulation or historical CSV required.
root<-'.'
files<-file.path(root,'tasks',sprintf('task_%03d.rds',1:150))
stopifnot(all(file.exists(files)));objects<-lapply(files,readRDS)
for(i in 1:150) {
 x<-objects[[i]];stopifnot(x$identity$task==i,x$identity$mc==500L,
 x$identity$dep==(i-1L)%/%50L+1L,nrow(x$table)==30L,nrow(x$grid)==400L)
 expected<-(((i-1L)%%50L)*10L+1L):(((i-1L)%%50L)*10L+10L)
 stopifnot(identical(as.integer(x$identity$reps),expected),
           identical(sort(unique(x$table$sim_id)),expected),
           identical(sort(unique(x$grid$sim_id)),expected))
}
tab<-do.call(rbind,lapply(objects,`[[`,'table'));grid<-do.call(rbind,lapply(objects,`[[`,'grid'))
stopifnot(nrow(tab)==4500L,nrow(grid)==60000L,
 !anyDuplicated(tab[c('dep_index','sim_id','design')]),
 !anyDuplicated(grid[c('dep_index','sim_id','eta','rho')]),
 all(table(tab$dep_index,tab$design)==500L))
methods<-c('LM_Crude','LM_Subclass','RI_Naive','Oracle_Normal','Oracle_Empirical')
stopifnot(all(is.finite(as.matrix(tab[methods]))),all(as.matrix(tab[methods])>=0),all(as.matrix(tab[methods])<=1))
out<-file.path(root,'results');dir.create(out,showWarnings=FALSE)
stopifnot(all(is.finite(grid$p_value)), all(grid$p_value >= 0 & grid$p_value <= 1),
          all(grid$design == 'P_2'))
saveRDS(list(table=tab,grid=grid),file.path(out,'null_rows.rds'))
summary<-do.call(rbind,lapply(split(tab,interaction(tab$dep_index,tab$design,drop=TRUE)),function(x) {
 data.frame(dep_index=x$dep_index[1],design=x$design[1],method=methods,
            rate=vapply(x[methods],function(z)mean(z<.05),numeric(1)),n=nrow(x))
}))
saveRDS(summary,file.path(out,'null_rejection_rates.rds'))
labels<-c('Independent','Moderate','Strong')
sens<-do.call(rbind,lapply(split(grid,interaction(grid$dep_index,grid$eta,grid$rho,drop=TRUE)),function(x) {
 stopifnot(nrow(x)==500L)
 data.frame(dep_index=x$dep_index[1],eta=x$eta[1],rho=x$rho[1],rate=mean(x$p_value<.05),n=nrow(x))
}))
sens$dependence<-factor(labels[sens$dep_index],levels=labels)
saveRDS(sens,file.path(out,'sensitivity_plot_data.rds'))
naive<-summary[summary$design=='P_2' & summary$method=='RI_Naive',]
naive$dependence<-factor(labels[naive$dep_index],levels=labels)
library(ggplot2)
p<-ggplot(sens,aes(rho,rate,colour=factor(eta),group=factor(eta)))+
 geom_hline(yintercept=.05,linetype='dashed',colour='black',linewidth=.7)+
 geom_hline(data=naive,aes(yintercept=rate),colour='red3',linetype='dotted',linewidth=.7)+
 geom_line(linewidth=.8)+geom_point(size=1.8)+
 scale_colour_manual(values=colorRampPalette(c('#2166AC','#74ADD1','#FEE090','#F46D43','#A50026'))(5),name=expression(eta))+
 scale_x_continuous(breaks=seq(.05,.40,.05),name=expression(rho))+
 scale_y_continuous(name='Rejection rate')+facet_wrap(~dependence,ncol=3)+theme_bw(base_size=12)+
 theme(strip.background=element_rect(fill='grey92',colour='grey60'),strip.text=element_text(face='bold'),panel.grid.minor=element_blank(),axis.text.x=element_text(angle=45,hjust=1),legend.position='right')
ggsave(file.path(out,'sensitivity_replication.pdf'),p,width=11,height=4.8)
ggsave(file.path(out,'sensitivity_replication.png'),p,width=11,height=4.8,dpi=180)
cat('PASS: 150 tasks, 1500 original matched datasets, 4500 oracle rows and 60000 grid rows.\n')
print(summary,row.names=FALSE)
