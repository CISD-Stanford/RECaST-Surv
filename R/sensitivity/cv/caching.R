# Selectively extracted definitions; no historical runner is sourced.
# See docs/PROVENANCE.csv for source identity and extraction scope.
O08_states <- function(chain) {
 stopifnot(is.matrix(chain),ncol(chain)==4L,all(is.finite(chain)))
 # Full IEEE double byte representation, not numeric rounding or a hash.
 keys<-vapply(seq_len(nrow(chain)),function(i)
  paste(sprintf('%02x',as.integer(writeBin(as.double(chain[i,]),raw(),size=8,endian='little'))),collapse=''),character(1))
 first<-which(!duplicated(keys));map<-match(keys,keys[first]);mult<-tabulate(map,nbins=length(first))
 stopifnot(sum(mult)==nrow(chain),all(keys==keys[first][map]))
 list(keys=keys,first=first,map=map,multiplicity=mult)
}

O08_gl_pair <- function(c0,a,n,lp) {
 c0<-unname(c0);a<-unname(a);lp<-unname(lp)
 za<-analytic(c0,a,'density',lp);zs<-analytic(c0,a,'survival',0)
 if(!is.null(za))return(list(za,zs))
 cuts<-logH_cuts(c0,a,-48,6);q<-gl_rule(n)
 half<-diff(cuts)/2;mid<-cuts[-length(cuts)]+half
 v<-sweep(outer(q$nodes,half,'*'),2,mid,'+')
 lw<-sweep(matrix(log(q$weights),nrow=n,ncol=length(half)),2,log(half),'+')
 x<-c0+v;h<-exp(x);common<-dcauchy(v,0,a,log=TRUE)+lw
 kd<-rep(-Inf,length(x));ok<-is.finite(x)&is.finite(h);kd[ok]<-x[ok]-h[ok]
 # Same positive sum in the same node order as the scalar implementation.
 ld<-lsum(kd+common);ls<-lsum(-h+common)
 massL<-pcauchy(-48,location=c0,scale=a)
 massR<-pcauchy(6,location=c0,scale=a,lower.tail=FALSE)
 if(massL>0)ls<-lsum(c(ls,log(massL)))
 left<-exp(-48)*massL
 list(finish(exp(ld),messages='OK_FIXED_SEGMENT_RULE_NO_INTRINSIC_ERROR_ESTIMATE',
  evaluations=length(v),logprefactor=lp,logvalue=ld,tail_bound=left+exp(6-exp(6))*massR),
  finish(exp(ls),messages='OK_FIXED_SEGMENT_RULE_NO_INTRINSIC_ERROR_ESTIMATE',
  evaluations=length(v),logvalue=ls,tail_bound=left+exp(-exp(6))*massR))
}

O08_pair <- function(c0,a,lp) {
 tm<-proc.time()[3];z64<-O08_gl_pair(c0,a,64,lp);s64<-unname(proc.time()[3]-tm)
 tm<-proc.time()[3];z128<-O08_gl_pair(c0,a,128,lp);s128<-unname(proc.time()[3]-tm)
 # Independent adaptive S/D calls retained; neither is reconstructed from GL.
 lapply(1:2,function(b) {
  comp<-c('density','survival')[b];pref<-if(b==1)lp else 0
  tm<-proc.time()[3]
  ad<-centered_adaptive(unname(c0),unname(a),comp,rel=1e-11,L=-48,U=6,logprefactor=unname(pref))
  sad<-unname(proc.time()[3]-tm);z<-z64[[b]];r<-z128[[b]]
  metrics<-c(log64=z$logvalue,log128=r$logvalue,log_adaptive=ad$logvalue,
   GL_refinement=abs(z$logvalue-r$logvalue),adaptive_difference=abs(z$logvalue-ad$logvalue),
   relative_tail=max(z$tail_bound/z$value_dimensionless,r$tail_bound/r$value_dimensionless,ad$tail_bound/ad$value_dimensionless),
   relative_adaptive_error=ad$absolute_error_estimate/ad$value_dimensionless,
   nodes64=z$evaluations,nodes128=r$evaluations,adaptive_subintervals=ad$evaluations,
   seconds64=s64/2,seconds128=s128/2,seconds_adaptive=sad)
  ok<-all(is.finite(metrics))&&metrics['GL_refinement']<=1e-8&&metrics['adaptive_difference']<=1e-8&&
   metrics['relative_tail']<=1e-10&&metrics['relative_adaptive_error']<=1e-8&&ad$messages%in%c('OK','ANALYTIC_ZERO_SCORE')
  list(ok=unname(ok),metrics=metrics,GL64=z,GL128=r,adaptive=ad)
 })
}

O08_score <- function(chain,time,score,ids,cache=FALSE,joint=FALSE,state=new.env(parent=emptyenv())) {
 stopifnot(is.matrix(chain),ncol(chain)==4L,nrow(chain)==801L,all(is.finite(chain)),
  length(time)==length(score),length(ids)==length(time),all(is.finite(time)&time>0),
  all(is.finite(score)),all(is.finite(exp(chain[,2:3]))&exp(chain[,2:3])>0))
 M<-nrow(chain);N<-length(time);states<-O08_states(chain)
 first<-if(cache)states$first else seq_len(M);map<-if(cache)states$map else seq_len(M)
 U<-length(first);inner<-array(NA_real_,c(U,N,2,3));checks<-vector('list',U*N*2);k<-0L
 for(i in seq_len(N))for(u in seq_len(U)) {
  m<-first[u];par<-unname(chain[m,]);c0<-exp(par[3])*log(time[i])+par[4]+par[1]*score[i]
  a<-exp(par[2])*abs(score[i]);lp<-par[3]-log(time[i])
  context<-list(patient_id=ids[i],patient_index=i,posterior_row=m,time=time[i],score=score[i],parameter=par,c=c0,a=a)
  state$context<-context
  z<-tryCatch(if(joint)O08_pair(c0,a,lp) else
   list(A08_one(c0,a,'density',lp),A08_one(c0,a,'survival',0)),error=function(e) {
    state$failure<-c(context,list(error=conditionMessage(e)));stop(e)
   })
  for(b in 1:2) {
   if(!isTRUE(z[[b]]$ok)) {
    state$failure<-c(context,list(component=b,result=z[[b]],raw=A08_failure_details(c0,a,c('density','survival')[b],if(b==1)lp else 0)))
    stop('Unchanged numerical gate failed; no fallback')
   }
   k<-k+1L;checks[[k]]<-c(patient_index=i,posterior_row=m,component=b,c=c0,a=a,z[[b]]$metrics)
   inner[u,i,b,]<-z[[b]]$metrics[c('log64','log128','log_adaptive')]
  }
 }
 # Expand in ORIGINAL chronology before the identical 801-row log-sum-exp.
 expanded<-lapply(1:3,function(q)array(inner[map,,,q],c(M,N,2),
  dimnames=list(NULL,as.character(ids),c('logf','logS'))))
 av<-lapply(expanded,function(x)t(vapply(seq_len(N),function(i)vapply(1:2,function(b)
  lsum(x[,i,b])-log(M),numeric(1)),numeric(2))))
 for(q in 1:3)dimnames(av[[q]])<-list(as.character(ids),c('logf','logS'))
 gap<-max(abs(c(av[[1]]-av[[2]],av[[1]]-av[[3]])))
 stopifnot(all(is.finite(unlist(av))),gap<=1e-4,all(av[[1]][,2]<=0),all(is.finite(exp(av[[1]]))),all(exp(av[[1]])>0))
 list(log_predictions=av[[1]],log_predictions128=av[[2]],log_predictions_adaptive=av[[3]],
  inner64=expanded[[1]],inner128=expanded[[2]],inner_adaptive=expanded[[3]],
  diagnostics=do.call(rbind,checks),state_map=states,evaluated_rows=first,
  posterior_predictive_max_log_gap=gap,all_adaptive_messages_OK=TRUE,
  conceptual_scalar_checks=2*M*N,evaluated_scalar_checks=2*U*N,
  GL_rule_calls=if(joint)2*U*N else 4*U*N,adaptive_calls=2*U*N)
}

