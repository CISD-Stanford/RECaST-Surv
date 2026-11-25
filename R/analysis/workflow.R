release_reference <- function(scenario,method) {
  if(scenario$family=='main') {
    r<-read.csv(file.path(release_root,'config/main/references.csv'),stringsAsFactors=FALSE)
    z<-r[r$scenario_id==scenario$id & r$method==method,,drop=FALSE]
    if(nrow(z)==1L)return(list(value=z$reference,id=z$reference_id))
    # No new population Cox projection is calculated for extension cells.
    if(method=='Cox-A')return(list(value=exp(scenario$spec$delta),id='conditional_PH_HR'))
    if(method=='Cox-U' && scenario$spec$delta==0)return(list(value=1,id='null_Cox_reference'))
  }
  if(scenario$family=='H1-H5') {
    r<-read.csv(file.path(release_root,'config/H1-H5/references.csv'),stringsAsFactors=FALSE)
    name<-switch(method,median='oracle_marginal_OE',mean='oracle_marginal_OE',
      'Cox-U'='treatment_only_population_Cox_HR','Cox-A'='conditional_adjusted_PH_Cox_HR')
    z<-r[r$scenario==scenario$spec$scenario & r$effect_id==scenario$spec$variant & r$reference_id==name,,drop=FALSE]
    if(nrow(z)==1L && z$applicability=='APPLICABLE')return(list(value=as.numeric(z$authoritative_value),id=name))
  }
  list(value=NA_real_,id='NO_APPROVED_REFERENCE_BOUND')
}

release_result_rows <- function(x,scenario,replication_id,smoke) {
  result<-list()
  for(m in names(x$calibration)) {
    z<-x$final$methods[[m]]
    s<-extract_method_splits(x$splits,m)
    ci<-construct_ci(z$O,z$E,s$O,s$E)
    ref<-release_reference(scenario,m)
    for(cal in c(TRUE,FALSE)) {
      result[[length(result)+1L]]<-data.frame(scenario_id=scenario$id,canonical_id=scenario$canonical_id,
        replication_id=replication_id,method=paste0('RECaST_',m,if(cal)'_calibrated' else '_uncalibrated'),
        point_estimator=m,estimate=z$HR,p=z$p_value,cutoff=if(cal)x$calibration[[m]]$threshold else .05,
        reject=if(cal)x$calibration[[m]]$reject else z$p_value<.05,
        O=z$O,E=z$E,CI_status=ci$status,CI_lower=ci$lower,CI_upper=ci$upper,CI_length=ci$length,
        reference=ref$value,reference_id=ref$id,smoke_only=smoke,stringsAsFactors=FALSE)
    }
  }
  do.call(rbind,result)
}

release_run_one <- function(scenario,state,replication_id,out,smoke=FALSE,generation_only=FALSE) {
  stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
  assign('.Random.seed',as.integer(state),.GlobalEnv)
  started<-proc.time();warnings<-character();status<-'FAILED';error<-NULL;x<-NULL;generated<-NULL
  tracker<-new.env(parent=emptyenv());tracker$stages<-list();tracker$fit_records<-list()
  ctx<-release_context(scenario,smoke)
  result<-tryCatch(withCallingHandlers({
    generated<-release_generate(scenario,ctx)
    saveRDS(generated,file.path(out,'inputs.rds'))
    if(generation_only) {status<-'GENERATED_ONLY';NULL} else {
    models<-ctx$engine$fit_external_models(generated$external)
    coefs<-if(scenario$family=='A-D')AD_source_coefficients(models,ctx$config) else lapply(models,coef)
    x<-ctx$functions$A03a_analyze_pair(generated$trial,coefs,ctx$config,ctx$engine,tracker)
    x$inputs<-generated;x$source_models<-models
    x$metadata<-list(config=ctx$config,delta=ctx$config$delta,scenario_id=scenario$id,replication_id=replication_id,
      rng_start=state,smoke_only=smoke)
    stopifnot(all(vapply(x$calibration,function(z)z$n_nonfinite==0L,logical(1))))
    for(st in c(x$splits,list(x$final))) {
      stopifnot(all(vapply(st$fit_diagnostics,function(z)z$fitting_status=='returned_finite_chain' &&
        z$conditional_survival_checks$finite && z$conditional_survival_checks$in_unit_interval &&
        z$conditional_survival_checks$nonincreasing,logical(1))))
      for(z in st$methods)stopifnot(is.finite(z$E),z$E>0,all(is.finite(unlist(z[c('O','HR','Z','p_value')]))))
    }
    saveRDS(x,file.path(out,'analysis.rds'))
    tab<-release_result_rows(x,scenario,replication_id,smoke)
    td<-generated$trial$data
    trial<-list(data=data.frame(time=td$time,event=td$status,treatment=td$trt,td[c('X1','X2','X3')]),
      trial_id=paste0(scenario$canonical_id,'/',replication_id),bank_id=scenario$canonical_id,
      replication_id=replication_id,data_hash=digest::digest(generated$trial,algo='sha256'),source_identity='local_release_new_run')
    fits<-setNames(lapply(c('Cox-U','Cox-A'),function(m)fit_once(trial,m)),c('Cox-U','Cox-A'))
    saveRDS(fits,file.path(out,'cox.rds'))
    for(m in names(fits)) {
      z<-fits[[m]];ref<-release_reference(scenario,m);valid<-isTRUE(z$inference_available)
      row<-tab[1,,drop=FALSE];row$method<-m;row$point_estimator<-m
      row$estimate<-if(valid)z$HR else NA_real_;row$p<-if(valid)z$p else NA_real_
      row$cutoff<-.05;row$reject<-if(valid)z$reject else NA
      row$O<-row$E<-NA_real_;row$CI_status<-if(valid)'nonempty' else 'unavailable'
      row$CI_lower<-if(valid)z$CI_lower else NA_real_;row$CI_upper<-if(valid)z$CI_upper else NA_real_
      row$CI_length<-row$CI_upper-row$CI_lower;row$reference<-ref$value;row$reference_id<-ref$id
      tab<-rbind(tab,row)
    }
    write.csv(tab,file.path(out,'trial_results.csv'),row.names=FALSE,na='NA')
    status<-'COMPLETED';tab
    }
  },warning=function(w) {
    warnings<<-c(warnings,conditionMessage(w));stop('Scientific warning retained without retry: ',conditionMessage(w))
  }),error=function(e){error<<-conditionMessage(e);NULL})
  # Do not discard a failed identity or silently omit its denominator.
  record<-list(scenario_id=scenario$id,canonical_id=scenario$canonical_id,replication_id=replication_id,
    status=status,error=error,warnings=warnings,retries=0L,rng_start=state,rng_end=.Random.seed,
    elapsed_seconds=unname((proc.time()-started)[3]),smoke_only=smoke,config=ctx$config)
  saveRDS(record,file.path(out,'status.rds'))
  if(status=='FAILED')saveRDS(list(partial_stages=tracker$stages,fit_records=tracker$fit_records),file.path(out,'partial.rds'))
  record
}
