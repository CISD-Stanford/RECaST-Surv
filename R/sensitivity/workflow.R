# These adapters reuse the approved statistical kernels on newly saved outputs.
release_focus_scope <- function(x) {
  c<-x$metadata$config
  stopifnot(!isTRUE(x$metadata$smoke_only),c$trial_N==300,c$r==2,c$family=='Same',
    x$metadata$delta %in% c(0,-.5),length(x$splits)==100L,identical(c$times_grid,1:60))
}
release_reweight <- function(x,rule=c('original','equal','normalized','cv5'),fixed_weights=NULL) {
  release_focus_scope(x);rule<-match.arg(rule)
  stages<-c(list(x$final),x$splits);norms<-vapply(x$source_coefs,function(b)sqrt(sum(b^2)),numeric(1))
  val<-array(NA_real_,c(101,2,5),dimnames=list(NULL,c('median','mean'),c('O','E','HR','T','p')))
  weights<-matrix(NA_real_,101,3)
  for(b in 1:101) {
    st<-stages[[b]]
    w<-switch(rule,original={a<-1/(st$source_scales+1e-6);a/sum(a)},equal=rep(1/3,3),
      normalized={d<-st$source_scales*norms;stopifnot(all(is.finite(d)&d>0));lw<--log(st$source_scales)-log(norms);a<-exp(lw-max(lw));a/sum(a)},cv5=fixed_weights)
    stopifnot(length(w)==3L,all(is.finite(w)),all(w>=0),abs(sum(w)-1)<1e-12);weights[b,]<-w
    for(m in c('median','mean')) {
      curve<-drop(st$methods[[m]]$source_cohort_survival%*%w)
      stopifnot(all(is.finite(curve)&curve>0&curve<=1))
      h<-approx(1:60,-log(curve),xout=st$obs_time,rule=2)$y
      O<-sum(st$obs_status);E<-sum(h);T<-(O-E)/sqrt(E)
      val[b,m,]<-c(O,E,O/E,T,pnorm(T))
    }
  }
  inference<-setNames(lapply(c('median','mean'),function(m) {
    z<-val[,m,];A07_calibrate(z[1,'O'],z[1,'E'],z[-1,'O'],z[-1,'E'],100L,z[1,'p'])
  }),c('median','mean'))
  list(rule=rule,weights=weights,statistics=val,inference=inference)
}
release_cv5 <- function(x) {
  release_focus_scope(x)
  e<-environment(release_cv5)
  for(p in c('integration.R','reference_scorer.R','caching.R','scorer.R','optimizer_original.R','optimizer.R'))
    sys.source(file.path(release_root,'R/sensitivity/cv',p),e)
  tr<-x$inputs$trial;ids<-which(tr$trt==0);folds<-sample(rep(1:5,length.out=length(ids)))
  tab<-matrix(NA_real_,length(ids),3);fits<-list()
  for(b in 1:5)for(j in 1:3) {
    train<-ids[folds!=b];held<-ids[folds==b]
    fit<-mcmc_survival_ph(tr$time[train],tr$status[train],
      drop(tr$X_mat[train,,drop=FALSE]%*%x$source_coefs[[j]]),
      init_par=x$metadata$config$init_par,steps=1000L,burnin=200L)
    stopifnot(nrow(fit$chain)==801L,all(is.finite(fit$chain)),length(O08_states(fit$chain)$first)>1)
    lp<-vapply(1:801,function(i)log_post_survival_ph(fit$chain[i,],tr$time[train],tr$status[train],
      drop(tr$X_mat[train,,drop=FALSE]%*%x$source_coefs[[j]])),numeric(1))
    stopifnot(all(is.finite(lp)))
    pred<-F08_score(fit$chain,tr$time[held],tr$status[held],
      drop(tr$X_mat[held,,drop=FALSE]%*%x$source_coefs[[j]]),held)
    tab[match(held,ids),j]<-pred$log_required
    fits[[paste(b,j,sep='_')]]<-list(fit=fit,score=pred)
  }
  stopifnot(all(is.finite(tab)))
  opt<-stack_optimize(tab);stopifnot(opt$status=='unique_numerical_solution')
  list(folds=data.frame(id=ids,fold=folds),fits=fits,optimizer=opt,
    recombination=release_reweight(x,'cv5',opt$weights))
}
release_calibration_sensitivity <- function(x) {
  release_focus_scope(x);stopifnot(x$metadata$delta==0)
  tr<-x$inputs$trial;controls<-which(tr$trt==0);n0<-length(controls);n1<-sum(tr$trt==1)
  engine<-release_context(release_scenario(x$metadata$scenario_id))$engine
  extra<-lapply(1:400,function(b) {
    train<-controls[sample.int(n0,floor(n0/3),replace=FALSE)]
    A03a_stage(tr,x$source_coefs,train,setdiff(controls,train),x$metadata$config,engine,
      stage_id=sprintf('additional_split_%03d',b))
  })
  bootstrap<-lapply(1:100,function(b) {
    train<-controls[sample.int(n0,n0,replace=TRUE)]
    evaluate<-controls[sample.int(n0,n1,replace=TRUE)]
    A03a_stage(tr,x$source_coefs,train,evaluate,x$metadata$config,engine,
      stage_id=sprintf('bootstrap_%03d',b))
  })
  bank<-c(x$splits,extra)
  infer<-function(stages,B) setNames(lapply(c('median','mean'),function(m) {
    z<-x$final$methods[[m]];s<-extract_method_splits(stages,m)
    A07_calibrate(z$O,z$E,s$O,s$E,B,z$p_value)
  }),c('median','mean'))
  list(split_bank=bank,bootstrap=bootstrap,B100=infer(bank[1:100],100L),
    B200=infer(bank[1:200],200L),B500=infer(bank,500L),bootstrap_B100=infer(bootstrap,100L))
}
