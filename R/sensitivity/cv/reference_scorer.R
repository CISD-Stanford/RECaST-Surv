# Selectively extracted definitions; no historical runner is sourced.
# See docs/PROVENANCE.csv for source identity and extraction scope.
A08_one <- function(c0,a,component,lp=0) {
 c0<-unname(c0);a<-unname(a);lp<-unname(lp)
 tm<-Sys.time();z64<-centered_gl(c0,a,component,64,logprefactor=lp)
 s64<-as.numeric(difftime(Sys.time(),tm,units='secs'))
 tm<-Sys.time();z128<-centered_gl(c0,a,component,128,logprefactor=lp)
 s128<-as.numeric(difftime(Sys.time(),tm,units='secs'))
 tm<-Sys.time();ad<-centered_adaptive(c0,a,component,rel=1e-11,L=-48,U=6,logprefactor=lp)
 sad<-as.numeric(difftime(Sys.time(),tm,units='secs'))
 metrics<-c(log64=z64$logvalue,log128=z128$logvalue,log_adaptive=ad$logvalue,
  GL_refinement=abs(z64$logvalue-z128$logvalue),
  adaptive_difference=abs(z64$logvalue-ad$logvalue),
  relative_tail=max(z64$tail_bound/z64$value_dimensionless,
   z128$tail_bound/z128$value_dimensionless,ad$tail_bound/ad$value_dimensionless),
  relative_adaptive_error=ad$absolute_error_estimate/ad$value_dimensionless,
  nodes64=z64$evaluations,nodes128=z128$evaluations,adaptive_subintervals=ad$evaluations,
  seconds64=s64,seconds128=s128,seconds_adaptive=sad)
 analytic_ok<-ad$messages=='ANALYTIC_ZERO_SCORE'
 ok<-all(is.finite(metrics))&&metrics['GL_refinement']<=1e-8&&
  metrics['adaptive_difference']<=1e-8&&metrics['relative_tail']<=1e-10&&
  metrics['relative_adaptive_error']<=1e-8&&(ad$messages=='OK'||analytic_ok)
 list(ok=unname(ok),metrics=metrics,GL64=z64,GL128=z128,adaptive=ad)
}

A08_failure_details <- function(c0,a,component,lp) {
 # Deterministic raw node/cut reconstruction only, NOT a scoring retry.
 tryCatch({
  cuts<-if(a>0)logH_cuts(c0,a) else numeric()
  nodes<-lapply(c(64,128),function(n) {
   if(a==0)return(list(analytic=TRUE))
   q<-gl_rule(n);half<-diff(cuts)/2;mid<-cuts[-length(cuts)]+half
   v<-sweep(outer(q$nodes,half,'*'),2,mid,'+')
   list(v=v,logH=c0+v,weights=q$weights,halfwidths=half,
    log_kernel=kernel_log(c0+v,component),log_cauchy_density=dcauchy(v,0,a,log=TRUE))
  });list(cuts=cuts,raw_nodes=nodes,logprefactor=lp)
 },error=function(e)list(raw_reconstruction_error=conditionMessage(e)))
}

A08_score_chain <- function(chain,time,score,ids,state=new.env(parent=emptyenv())) {
 stopifnot(is.matrix(chain),ncol(chain)==4L,all(is.finite(chain)),nrow(chain)==801L,
  length(time)==length(score),length(ids)==length(time),all(is.finite(time)&time>0),
  all(is.finite(score)),all(is.finite(exp(chain[,2:3]))&exp(chain[,2:3])>0))
 M<-nrow(chain);N<-length(time);state$completed_checks<-0L;state$failure<-NULL
 fields<-c('patient_index','posterior_row','component','c','a','log64','log128','log_adaptive',
  'GL_refinement','adaptive_difference','relative_tail','relative_adaptive_error',
  'nodes64','nodes128','adaptive_subintervals','seconds64','seconds128','seconds_adaptive')
 state$checks<-matrix(NA_real_,2*M*N,length(fields),dimnames=list(NULL,fields))
 state$inner64<-state$inner128<-state$inner_adaptive<-array(NA_real_,c(M,N,2),
  dimnames=list(NULL,as.character(ids),c('logf','logS')))
 for(i in seq_len(N))for(m in seq_len(M))for(b in 1:2) {
  par<-unname(chain[m,]);c0<-exp(par[3])*log(time[i])+par[4]+par[1]*score[i]
  a<-exp(par[2])*abs(score[i]);comp<-c('density','survival')[b]
  lp<-if(b==1)par[3]-log(time[i]) else 0
  context<-list(patient_id=ids[i],patient_index=i,posterior_row=m,time=time[i],score=score[i],
   parameter=par,c=c0,a=a,component=comp,logprefactor=lp)
  z<-tryCatch(A08_one(c0,a,comp,lp),error=function(e)list(ok=FALSE,error=conditionMessage(e)))
  if(!isTRUE(z$ok)) {
   state$failure<-c(context,list(result=z,raw=A08_failure_details(c0,a,comp,lp)))
   stop('Approved numerical scoring gate failed at patient ',ids[i],', row ',m,', ',comp)
  }
  k<-state$completed_checks+1L
  state$checks[k,]<-c(i,m,b,c0,a,z$metrics)
  state$inner64[m,i,b]<-z$metrics['log64'];state$inner128[m,i,b]<-z$metrics['log128']
  state$inner_adaptive[m,i,b]<-z$metrics['log_adaptive'];state$completed_checks<-k
 }
 average<-function(x)t(vapply(seq_len(N),function(i)vapply(1:2,function(b)
  lsum(x[,i,b])-log(M),numeric(1)),numeric(2)))
 a64<-average(state$inner64);a128<-average(state$inner128);ad<-average(state$inner_adaptive)
 dimnames(a64)<-dimnames(a128)<-dimnames(ad)<-list(as.character(ids),c('logf','logS'))
 gap<-max(abs(c(a64-a128,a64-ad)))
 if(any(!is.finite(c(a64,a128,ad)))||gap>1e-4||any(a64[,'logS']>0)||
    any(!is.finite(exp(a64)))||any(exp(a64)<=0)) {
  state$failure<-list(reason='Unusable posterior predictive average',log64=a64,log128=a128,
   log_adaptive=ad,average_log_gap=gap);stop('Posterior predictive log-score gate failed')
 }
 list(log_predictions=a64,log_predictions128=a128,log_predictions_adaptive=ad,
  f=exp(a64[,'logf']),S=exp(a64[,'logS']),M=M,ids=ids,
  posterior_predictive_max_log_gap=gap,all_adaptive_messages_OK=TRUE,
  diagnostics=state$checks,inner64=state$inner64,inner128=state$inner128,
  inner_adaptive=state$inner_adaptive)
}

