release_summarize <- function(input,out) {
  stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
  manifest<-readRDS(file.path(input,'run_manifest.rds'))
  tables<-list();trial_tables<-list()
  for(s in manifest$memberships) {
    files<-file.path(input,s$canonical_id,sprintf('rep_%04d',seq_len(sum(manifest$plan$canonical_id==s$canonical_id))),'trial_results.csv')
    available<-files[file.exists(files)]
    dat<-if(length(available))do.call(rbind,lapply(available,read.csv,stringsAsFactors=FALSE)) else NULL
    methods<-c(paste0('RECaST_',if(s$family=='A-D')'median' else rep(c('median','mean'),each=2),
      if(s$family=='A-D')c('_calibrated','_uncalibrated') else rep(c('_calibrated','_uncalibrated'),2)), 'Cox-U','Cox-A')
    for(m in methods) {
      z<-if(is.null(dat))NULL else dat[dat$method==m,,drop=FALSE]
      planned<-length(files);point_n<-if(is.null(z))0L else sum(is.finite(z$estimate))
      test_n<-if(is.null(z))0L else sum(!is.na(z$reject));k<-if(is.null(z))0L else sum(z$reject,na.rm=TRUE)
      rate<-if(test_n)k/test_n else NA_real_
      avg<-function(v)if(length(v))mean(v)else NA_real_
      mc<-function(v)if(length(v)>1L)sd(v)/sqrt(length(v))else NA_real_
      estimates<-if(is.null(z))numeric()else z$estimate[is.finite(z$estimate)]
      pm<-sub('^RECaST_','',sub('_(un)?calibrated$','',m))
      ref<-release_reference(s,pm)
      eligible<-if(is.null(z))logical()else z$CI_status %in% c('nonempty','empty')
      eligible_n<-sum(eligible);nonempty<-if(is.null(z))logical()else z$CI_status=='nonempty'
      contain<-if(!is.finite(ref$value)||is.null(z))logical()else ifelse(z$CI_status[eligible]=='empty',FALSE,
        z$CI_lower[eligible]<=ref$value & ref$value<=z$CI_upper[eligible])
      low<-if(length(contain))ifelse(z$CI_status[eligible]=='empty',FALSE,z$CI_upper[eligible]<ref$value)else logical()
      up<-if(length(contain))ifelse(z$CI_status[eligible]=='empty',FALSE,z$CI_lower[eligible]>ref$value)else logical()
      widths<-if(is.null(z))numeric()else z$CI_length[nonempty & is.finite(z$CI_length)]
      errors<-if(is.finite(ref$value))estimates-ref$value else numeric()
      tables[[length(tables)+1L]]<-data.frame(scenario_id=s$id,canonical_id=s$canonical_id,method=m,
        planned=planned,represented=if(is.null(z))0L else nrow(z),point_n=point_n,test_n=test_n,
        unknown_test_n=planned-test_n,rejections=k,rejection_rate=rate,
        rejection_MCSE=if(test_n)sqrt(rate*(1-rate)/test_n)else NA_real_,
        planned_rate_lower=k/planned,planned_rate_upper=(k+planned-test_n)/planned,
        mean_estimate=avg(estimates),empirical_SD=if(point_n>1L)sd(estimates)else NA_real_,
        mean_MCSE=mc(estimates),reference=ref$value,reference_id=ref$id,
        mean_error=avg(errors),MSE=avg(errors^2),error_MCSE=mc(errors),MSE_MCSE=mc(errors^2),
        interval_eligible_n=eligible_n,interval_empty_n=if(is.null(z))0L else sum(z$CI_status=='empty'),
        interval_unavailable_n=planned-eligible_n,coverage_n=length(contain),
        coverage=avg(contain),coverage_MCSE=if(length(contain))sqrt(mean(contain)*(1-mean(contain))/length(contain))else NA_real_,
        lower_tail_noncoverage=avg(low),upper_tail_noncoverage=avg(up),
        length_n=length(widths),mean_length=avg(widths),length_MCSE=mc(widths),smoke_only=isTRUE(manifest$smoke_only))
    }
    if(!is.null(dat)) {dat$scenario_id<-s$id;trial_tables[[length(trial_tables)+1L]]<-dat}
  }
  tab<-do.call(rbind,tables);write.csv(tab,file.path(out,'summary.csv'),row.names=FALSE,na='NA')
  if(length(trial_tables))write.csv(do.call(rbind,trial_tables),file.path(out,'trial_results.csv'),row.names=FALSE,na='NA')
  writeLines(c('New-run summaries only; no historical outcomes bundled.',
    'Rates use available tests; bounds retain all planned trials. Empty intervals are noncoverage; length uses nonempty intervals.',
    'Reference-NA quantities remain NA. Reused memberships are not independent trials.',
    if(isTRUE(manifest$smoke_only))'SMOKE ONLY: reduced computation; no paper inference.' else 'Production configuration used; numerical convergence is not established.'),file.path(out,'SUMMARY_NOTES.txt'))
  grDevices::pdf(file.path(out,'rejection_rates.pdf'),width=10,height=6)
  for(id in unique(tab$scenario_id)) {
    z<-tab[tab$scenario_id==id,]
    if(all(!is.finite(z$rejection_rate)))next
    par(mar=c(11,4,4,1))
    labels<-sub('^RECaST_','',z$method)
    labels<-gsub('_',', ',labels,fixed=TRUE)
    barplot(z$rejection_rate,names.arg=labels,cex.names=.8,las=2,ylim=c(0,1),ylab='Available-test rejection fraction',
      main=paste(if(isTRUE(manifest$smoke_only))'SMOKE ONLY -' else 'New simulation -',id))
    abline(h=.05,lty=2,col='red')
  }
  grDevices::dev.off();invisible(tab)
}
