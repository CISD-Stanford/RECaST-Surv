# Selectively extracted definitions; no historical runner is sourced.
# See docs/PROVENANCE.csv for source identity and extraction scope.
F08_rule <- local({cache<-new.env(parent=emptyenv());function(n) {
 key<-as.character(n)
 if(!exists(key,cache,inherits=FALSE)) {
  q<-gl_rule(n);assign(key,list(nodes=q$nodes,logweights=log(q$weights),n=n),cache)
 }
 get(key,cache,inherits=FALSE)
}})

F08_one <- function(c0,a,component,lp=0) {
 c0<-unname(c0);a<-unname(a);lp<-unname(lp)
 az<-analytic(c0,a,component,lp)
 if(is.null(az)) {
  cuts<-logH_cuts(c0,a,-48,6);tails<-tail_terms(c0,a,component,-48,6)
  half<-diff(cuts)/2;mid<-cuts[-length(cuts)]+half;lh<-log(half)
  gl<-function(n) {
   q<-F08_rule(n);v<-sweep(outer(q$nodes,half,'*'),2,mid,'+')
   lw<-sweep(matrix(q$logweights,nrow=n,ncol=length(half)),2,lh,'+')
   # Preserve scalar reference summation order, including prefactor per row.
   lv<-lsum(c(kernel_log(c0+v,component)+dcauchy(v,0,a,log=TRUE)+lw))
   if(tails$add>0)lv<-lsum(c(lv,log(tails$add)))
   finish(exp(lv),messages='OK_FIXED_SEGMENT_RULE_NO_INTRINSIC_ERROR_ESTIMATE',
    evaluations=length(v),logprefactor=lp,logvalue=lv,tail_bound=tails$bound)
  }
  z64<-gl(64);z128<-gl(128)
  # Immutable independent adaptive implementation; no changed check frequency.
  ad<-centered_adaptive(c0,a,component,rel=1e-11,L=-48,U=6,logprefactor=lp)
 } else z64<-z128<-ad<-az
 metrics<-c(log64=z64$logvalue,log128=z128$logvalue,log_adaptive=ad$logvalue,
  GL_refinement=abs(z64$logvalue-z128$logvalue),adaptive_difference=abs(z64$logvalue-ad$logvalue),
  relative_tail=max(z64$tail_bound/z64$value_dimensionless,z128$tail_bound/z128$value_dimensionless,ad$tail_bound/ad$value_dimensionless),
  relative_adaptive_error=ad$absolute_error_estimate/ad$value_dimensionless,
  nodes64=z64$evaluations,nodes128=z128$evaluations,adaptive_subintervals=ad$evaluations)
 ok<-all(is.finite(metrics))&&metrics['GL_refinement']<=1e-8&&metrics['adaptive_difference']<=1e-8&&
  metrics['relative_tail']<=1e-10&&metrics['relative_adaptive_error']<=1e-8&&ad$messages%in%c('OK','ANALYTIC_ZERO_SCORE')
 list(ok=unname(ok),metrics=metrics,adaptive_message=ad$messages,
  failure_values=if(!isTRUE(ok))list(GL64=z64,GL128=z128,adaptive=ad) else NULL)
}

F08_score <- function(chain,time,status,score,ids,state=new.env(parent=emptyenv())) {
 stopifnot(is.matrix(chain),nrow(chain)==801L,ncol(chain)==4L,all(is.finite(chain)),
  length(time)==length(status),length(time)==length(score),length(ids)==length(time),
  all(is.finite(time)&time>0),all(status%in%c(0,1)),all(is.finite(score)),
  all(is.finite(exp(chain[,2:3]))&exp(chain[,2:3])>0))
 sm<-O08_states(chain);ch<-unname(chain[sm$first,,drop=FALSE]);U<-nrow(ch);N<-length(time)
 nu<-exp(ch[,3]);tau<-exp(ch[,2]);lt<-log(time);component<-ifelse(status==1,'density','survival')
 fields<-c('patient_index','unique_state','posterior_row','component','c','a','logprefactor',
  'log64','log128','log_adaptive','GL_refinement','adaptive_difference','relative_tail',
  'relative_adaptive_error','nodes64','nodes128','adaptive_subintervals')
 checks<-matrix(NA_real_,U*N,length(fields),dimnames=list(NULL,fields))
 inner<-array(NA_real_,c(U,N,3),dimnames=list(NULL,as.character(ids),c('GL64','GL128','adaptive')))
 messages<-matrix('',U,N);k<-0L
 for(i in seq_len(N)) {
  cv<-nu*lt[i]+ch[,4]+ch[,1]*score[i];av<-tau*abs(score[i])
  pref<-if(status[i]==1)ch[,3]-lt[i] else rep(0,U)
  for(u in seq_len(U)) {
   context<-list(patient_id=ids[i],patient_index=i,posterior_row=sm$first[u],unique_state=u,
    time=time[i],status=status[i],score=score[i],parameter=ch[u,],c=cv[u],a=av[u],
    component=component[i],logprefactor=pref[u])
   state$context<-context
   z<-tryCatch(F08_one(cv[u],av[u],component[i],pref[u]),error=function(e)list(ok=FALSE,error=conditionMessage(e)))
   if(!isTRUE(z$ok)) {
    state$failure<-c(context,list(result=z,raw=A08_failure_details(cv[u],av[u],component[i],pref[u])))
    stop('Required-component numerical gate failed; no retry/fallback')
   }
   k<-k+1L;checks[k,]<-c(i,u,sm$first[u],status[i],cv[u],av[u],pref[u],z$metrics)
   inner[u,i,]<-z$metrics[1:3];messages[u,i]<-z$adaptive_message
  }
 }
 averages<-vapply(seq_len(N),function(i)vapply(1:3,function(q)
  lsum(inner[sm$map,i,q])-log(801),numeric(1)),numeric(3))
 averages<-t(averages);dimnames(averages)<-list(as.character(ids),c('GL64','GL128','adaptive'))
 gap<-max(abs(c(averages[,1]-averages[,2],averages[,1]-averages[,3])))
 if(any(!is.finite(averages))||gap>1e-4||any(averages[status==0,1]>0)||
    any(!is.finite(exp(averages[,1])))||any(exp(averages[,1])<=0)) {
  state$failure<-list(reason='Unusable posterior average',averages=averages,inner=inner,state_map=sm)
  stop('Posterior-average gate failed; no repair')
 }
 weighted<-vapply(seq_len(N),function(i)lsum(log(sm$multiplicity)+inner[,i,1])-log(801),numeric(1))
 stopifnot(max(abs(weighted-averages[,1]))<=1e-12)
 list(ids=ids,time=time,status=status,source_score=score,requested_component=component,
  unused_component_status=rep('not requested',N),log_required=averages[,1],
  log_required_all_checks=averages,unique_log_contributions=inner,state_map=sm,
  checks=checks,adaptive_messages=messages,M=801L,posterior_max_log_gap=gap,
  conceptual_required_checks=801*N,evaluated_required_checks=U*N,
  GL_rule_calls=2*U*N,adaptive_calls=U*N,weighted_average_error=max(abs(weighted-averages[,1])))
}

