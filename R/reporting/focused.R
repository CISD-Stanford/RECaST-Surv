# Saved-new-run reporting formulas from accepted A06 v002 and A08c.
# The original complete 1000-per-effect designs are required; no missing-case
# grouping or new within-trial source uncertainty propagation is introduced.
release_report_focused <- function(input,out,kind=c('full-vs-split','coefficients')) {
  kind<-match.arg(kind);stopifnot(!dir.exists(out))
  get_trials<-function(id) {
    files<-file.path(input,id,sprintf('rep_%04d',1:1000),'analysis.rds')
    stopifnot(all(file.exists(files)))
    lapply(files,function(p){x<-readRDS(p);release_focus_scope(x);x})
  }
  null<-get_trials('Same_N300_r2_D0');dir.create(out,recursive=TRUE)
  write_new<-function(x,name)write.csv(x,file.path(out,name),row.names=FALSE,na='NA')
  if(kind=='full-vs-split') {
    tr<-sr<-gr<-list();grid<-seq(-5,5,by=.25)
    for(i in 1:1000)for(m in c('median','mean')) {
      x<-null[[i]];z<-x$final$methods[[m]];s<-extract_method_splits(x$splits,m)
      T<-(s$O-s$E)/sqrt(s$E);full<-(z$O-z$E)/sqrt(z$E)
      stopifnot(all(is.finite(c(T,full))))
      less<-sum(T<full);ties<-sum(T==full)
      tr[[length(tr)+1L]]<-data.frame(replication_id=i,method=m,T_full=full,
        raw_threshold=x$calibration[[m]]$raw_threshold,threshold=x$calibration[[m]]$threshold,
        rank_mid=(less+.5*ties)/100)
      sr[[length(sr)+1L]]<-data.frame(replication_id=i,method=m,split=1:100,T_split=T)
      ff<-as.numeric(full<=grid);fs<-vapply(grid,function(t)mean(T<=t),numeric(1))
      gr[[length(gr)+1L]]<-data.frame(replication_id=i,method=m,grid=grid,full_indicator=ff,within_trial_split_cdf=fs,paired_difference=ff-fs)
    }
    tr<-do.call(rbind,tr);sr<-do.call(rbind,sr);gr<-do.call(rbind,gr)
    moments<-quants<-cdf<-list()
    for(m in c('median','mean'))for(dist in c('full','split_equal_trial')) {
      v<-if(dist=='full')tr$T_full[tr$method==m]else sr$T_split[sr$method==m]
      moments[[length(moments)+1L]]<-data.frame(method=m,distribution=dist,n_trials=1000,
        scores_per_trial=if(dist=='full')1 else 100,score_count=length(v),mean=mean(v),empirical_sd=sqrt(mean((v-mean(v))^2)))
      probs<-c(.025,.05,.5,.95,.975)
      quants[[length(quants)+1L]]<-data.frame(method=m,distribution=dist,probability=probs,quantile=unname(quantile(v,probs,type=7)))
    }
    for(m in c('median','mean'))for(g in grid) {
      z<-gr[gr$method==m&gr$grid==g,]
      cdf[[length(cdf)+1L]]<-data.frame(method=m,grid=g,n_trials=nrow(z),full_cdf=mean(z$full_indicator),
        split_cdf=mean(z$within_trial_split_cdf),paired_difference=mean(z$paired_difference),
        paired_MCSE=sd(z$paired_difference)/sqrt(nrow(z)),split_cdf_cluster_MCSE=sd(z$within_trial_split_cdf)/sqrt(nrow(z)))
    }
    write_new(tr,'trial_ranks_and_thresholds.csv');write_new(do.call(rbind,moments),'score_moments.csv')
    write_new(do.call(rbind,quants),'score_quantiles.csv');write_new(do.call(rbind,cdf),'score_cdfs.csv')
    return(invisible(NULL))
  }
  alternative<-get_trials('Same_N300_r2_Dm05');sources<-inference<-list()
  for(setting in c('null','alternative')) {
    bank<-if(setting=='null')null else alternative
    for(i in 1:1000) {
      x<-bank[[i]];t<-x$inputs$trial;ids<-t$data$id[t$data$trt==0]
      for(j in 1:3) {
        b<-x$source_coefs[[j]];theta<-x$metadata$config$external_beta[j,];d<-b-theta
        lp_difference<-drop(t$X_mat[ids,,drop=FALSE]%*%b)-drop(t$X_mat[ids,,drop=FALSE]%*%theta)
        U<-sqrt(mean(lp_difference^2));D<-sqrt(sum(d^2))
        sources[[length(sources)+1L]]<-data.frame(setting=setting,replication_id=i,source=j,
          d_X1=d[1],d_X2=d[2],d_X3=d[3],D=D,U=U,weight=x$final$source_weights[j],row.names=NULL)
      }
      for(m in c('median','mean')) {
        z<-x$final$methods[[m]];ref<-if(setting=='null')1 else .64080597031188
        inference[[length(inference)+1L]]<-data.frame(setting=setting,replication_id=i,method=m,
          absolute_deviation=abs(z$HR-ref),squared_deviation=(z$HR-ref)^2,reject=x$calibration[[m]]$reject)
      }
    }
  }
  s<-do.call(rbind,sources);inf<-do.call(rbind,inference);groups<-coeff<-corr<-quartiles<-list()
  for(setting in c('null','alternative')) {
    u<-vapply(1:1000,function(i)mean(s$U[s$setting==setting&s$replication_id==i]),numeric(1))
    group<-integer(1000);ord<-order(u,sprintf('%04d',1:1000),method='radix');group[ord]<-rep(1:4,each=250)
    groups[[length(groups)+1L]]<-data.frame(setting=setting,replication_id=1:1000,Ubar=u,group=group)
    for(j in 1:3) {
      z<-s[s$setting==setting&s$source==j,]
      for(k in c('X1','X2','X3')) {v<-z[[paste0('d_',k)]]
        coeff[[length(coeff)+1L]]<-data.frame(setting=setting,source=j,coordinate=k,n=1000,
          mean_signed_error=mean(v),empirical_SD=sd(v),RMSE=sqrt(mean(v^2)))}
      corr[[length(corr)+1L]]<-data.frame(setting=setting,source=j,method=NA_character_,
        association='source_U_vs_original_weight',rho=cor(z$U,z$weight,method='spearman'))
    }
    for(m in c('median','mean')) {
      z<-inf[inf$setting==setting&inf$method==m,];z$Ubar<-u[z$replication_id];z$group<-group[z$replication_id]
      corr[[length(corr)+1L]]<-data.frame(setting=setting,source=NA_integer_,method=m,
        association='Ubar_vs_absolute_benchmark_deviation',rho=cor(z$Ubar,z$absolute_deviation,method='spearman'))
      for(g in 1:4) {
        q<-z[z$group==g,];stopifnot(nrow(q)==250L)
        quartiles[[length(quartiles)+1L]]<-data.frame(setting=setting,method=m,group=g,n=250,
          mean_absolute_deviation=mean(q$absolute_deviation),mean_squared_deviation=mean(q$squared_deviation),rejection_rate=mean(q$reject))
      }
    }
  }
  write_new(s,'source_discrepancies.csv');write_new(do.call(rbind,groups),'Ubar_groups.csv')
  write_new(do.call(rbind,coeff),'coefficient_summaries.csv');write_new(do.call(rbind,corr),'associations.csv')
  write_new(do.call(rbind,quartiles),'group_summaries.csv')
}
