release_rmap <- function(inputs,model_file,jags_seed) {
  required<-jsonlite::read_json(file.path(release_root,'external/RMAP_MODEL_REQUIRED.json'))
  if(!file.exists(model_file))stop('Authorized exact JAGS model file required; see external/RMAP_MODEL_REQUIRED.json')
  stopifnot(digest::digest(file=model_file,algo='sha256')==required$sha256)
  for(p in c('rjags','coda'))if(!requireNamespace(p,quietly=TRUE))stop('Missing optional dependency: ',p)
  e<-new.env(parent=environment(release_rmap))
  sys.source(file.path(release_root,'R/analysis/rmap_interface.R'),e)
  e$jags.model<-function(...) {
    a<-list(...);a$file<-normalizePath(model_file)
    a$inits<-list(.RNG.name='base::Mersenne-Twister',.RNG.seed=as.integer(jags_seed))
    do.call(rjags::jags.model,a)
  }
  e$update<-stats::update;e$coda.samples<-rjags::coda.samples
  tr<-inputs$trial$data
  groups<-c(lapply(inputs$external,`[[`,'data'),list(tr[tr$trt==0,,drop=FALSE],tr[tr$trt==1,,drop=FALSE]))
  z<-lapply(groups,e$convert_ind_data_to_poi_data)
  mydata<-list(r=vapply(z,`[[`,numeric(1),'r'),e=vapply(z,`[[`,numeric(1),'e'))
  stopifnot(length(mydata$r)==5L,all(is.finite(mydata$e)&mydata$e>0))
  warnings<-character()
  samples<-withCallingHandlers(e$rMAP_surv_fit(mydata,parameters=e$params,niter=2000),warning=function(w) {
    warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')
  })
  allowed<-'Failed to set trace monitor for log_lambda.pred\nVariable log_lambda.pred not found'
  if(any(trimws(warnings)!=allowed))stop('Unexpected R-MAP warning; no retry: ',paste(warnings,collapse=' | '))
  m<-as.matrix(samples);stopifnot(nrow(m)==400L,all(is.finite(m)),colnames(m)[1]=='beta1')
  list(score=sum(m[,1]<0)/length(m[,1]),estimate=mean(exp(m[,1])),reject_RMAP=mean(m[,1]<0)>.95,
    samples=samples,aggregates=mydata,warnings=warnings,jags_seed=jags_seed,
    convergence_assessed=FALSE)
}
release_crmap <- function(null_scores,evaluation_scores) {
  stopifnot(is.numeric(null_scores),length(null_scores)==1000L,
    all(is.finite(null_scores)&null_scores>=0&null_scores<=1),
    all(is.finite(evaluation_scores)&evaluation_scores>=0&evaluation_scores<=1))
  cutoff<-sort(null_scores)[950L]
  list(cutoff=cutoff,rank=950L,null_n=1000L,operator='>=',
    null_selection_rejections=null_scores>=cutoff,evaluation_rejections=evaluation_scores>=cutoff,
    qualification='Scenario-level calibration; selection-sample null frequency is not independent validation')
}
