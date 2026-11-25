# Selectively extracted definitions; no historical runner is sourced.
# See docs/PROVENANCE.csv for source identity and extraction scope.
lsum <- function(x) {b<-max(x);if(b==-Inf)return(-Inf);b+log(sum(exp(x-b)))}

features <- c(-48,-40,-32,-24,-16,-12,-8,-4,-2,-1,0,1,2,3,4,5,6)

gl_cache <- new.env(parent=emptyenv())

gl_rule <- function(n) {
 key<-as.character(n)
 if(!exists(key,gl_cache,inherits=FALSE))assign(key,statmod::gauss.quad(n,'legendre'),gl_cache)
 get(key,gl_cache,inherits=FALSE)
}

kernel_log <- function(x,component) {
 h<-exp(x)
 if(component=='survival') -h else {
  z<-rep(-Inf,length(x));ok<-is.finite(h)&is.finite(x);z[ok]<-x[ok]-h[ok];z
 }
}

finish <- function(value,error=NA_real_,messages='OK',evaluations=NA_integer_,
                   logprefactor=0,logvalue=NULL,tail_bound=0) {
 if(is.null(logvalue))logvalue<-log(value)
 list(logvalue=logvalue+logprefactor,value_dimensionless=value,
  absolute_error_estimate=error,tail_bound=tail_bound,messages=messages,
  evaluations=evaluations,usable=is.finite(logvalue))
}

analytic <- function(c,a,component,logprefactor) {
 if(length(c)!=1||length(a)!=1||!all(is.finite(c(c,a,logprefactor)))||a<0||
    !component%in%c('density','survival'))stop('Invalid experimental integral input; no floor or repair')
 if(a!=0)return(NULL)
 logv<-kernel_log(c,component)
 finish(exp(logv),0,'ANALYTIC_ZERO_SCORE',0,logprefactor,logv)
}

uniform_gl <- function(c,a,component,n,logprefactor=0) {
 z<-analytic(c,a,component,logprefactor);if(!is.null(z))return(z)
 q<-gl_rule(n);u<-pi*q$nodes/2
 lv<-lsum(kernel_log(c+a*tan(u),component)+log(q$weights/2))
 finish(exp(lv),messages='OK_FIXED_RULE_NO_ERROR_ESTIMATE',evaluations=n,
  logprefactor=logprefactor,logvalue=lv)
}

segmented_u <- function(c,a,component,rel=1e-10,logprefactor=0) {
 z<-analytic(c,a,component,logprefactor);if(!is.null(z))return(z)
 cuts<-sort(unique(c(-pi/2,0,atan((features-c)/a),pi/2)))
 warnings<-character()
 vals<-lapply(seq_len(length(cuts)-1),function(i)withCallingHandlers(
  integrate(function(u)exp(kernel_log(c+a*tan(u),component))/pi,cuts[i],cuts[i+1],
   rel.tol=rel,abs.tol=0,subdivisions=2000L,stop.on.error=FALSE),
  warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}))
 value<-sum(vapply(vals,function(z)z$value,numeric(1)))
 err<-sum(vapply(vals,function(z)z$abs.error,numeric(1)))
 msgs<-paste(unique(c(vapply(vals,function(z)z$message,character(1)),warnings)),collapse=';')
 finish(value,err,msgs,sum(vapply(vals,function(z)z$subdivisions,numeric(1))),logprefactor)
}

logH_cuts <- function(c,a,L=-48,U=6) {
 lo<-L-c;hi<-U-c
 # Centering at c keeps the Cauchy peak at v=0 and avoids subtraction of
 # two nearly equal floating-point x and c values inside its density.
 exponent<-ceiling(log10(max(abs(c(lo,hi)))/a))
 if(!is.finite(exponent)||exponent>300)stop('Unsupported geometric exponent range; no silent cut truncation')
 geo<-if(exponent>=0)a*10^seq.int(0,exponent) else numeric()
 candidates<-c(lo,hi,0,features-c,geo,-geo)
 sort(unique(candidates[is.finite(candidates)&candidates>=lo&candidates<=hi]))
}

tail_terms <- function(c,a,component,L,U) {
 massL<-pcauchy(L,location=c,scale=a)
 massR<-pcauchy(U,location=c,scale=a,lower.tail=FALSE)
 # S left approximation F(L), error <= exp(L)*F(L).
 # Density dimensionless D left <= exp(L)*F(L).
 left_error<-exp(L)*massL
 right_error<-if(component=='survival')exp(-exp(U))*massR else exp(U-exp(U))*massR
 list(add=if(component=='survival')massL else 0,bound=left_error+right_error)
}

centered_adaptive <- function(c,a,component,rel=1e-10,L=-48,U=6,logprefactor=0) {
 z<-analytic(c,a,component,logprefactor);if(!is.null(z))return(z)
 cuts<-logH_cuts(c,a,L,U);tails<-tail_terms(c,a,component,L,U)
 warnings<-character()
 vals<-lapply(seq_len(length(cuts)-1),function(i)withCallingHandlers(
  integrate(function(v)exp(kernel_log(c+v,component)+dcauchy(v,0,a,log=TRUE)),
   cuts[i],cuts[i+1],rel.tol=rel,abs.tol=0,subdivisions=2000L,stop.on.error=FALSE),
  warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}))
 value<-tails$add+sum(vapply(vals,function(z)z$value,numeric(1)))
 error<-tails$bound+sum(vapply(vals,function(z)z$abs.error,numeric(1)))
 msgs<-paste(unique(c(vapply(vals,function(z)z$message,character(1)),warnings)),collapse=';')
 finish(value,error,msgs,sum(vapply(vals,function(z)z$subdivisions,numeric(1))),
  logprefactor,tail_bound=tails$bound)
}

centered_gl <- function(c,a,component,n=64,L=-48,U=6,logprefactor=0) {
 z<-analytic(c,a,component,logprefactor);if(!is.null(z))return(z)
 cuts<-logH_cuts(c,a,L,U);tails<-tail_terms(c,a,component,L,U);q<-gl_rule(n)
 half<-diff(cuts)/2;mid<-cuts[-length(cuts)]+half
 v<-sweep(outer(q$nodes,half,'*'),2,mid,'+')
 lw<-sweep(matrix(log(q$weights),nrow=n,ncol=length(half)),2,log(half),'+')
 lv<-lsum(c(kernel_log(c+v,component)+dcauchy(v,0,a,log=TRUE)+lw))
 if(tails$add>0)lv<-lsum(c(lv,log(tails$add)))
 finish(exp(lv),messages='OK_FIXED_SEGMENT_RULE_NO_INTRINSIC_ERROR_ESTIMATE',
  evaluations=length(v),logprefactor=logprefactor,logvalue=lv,tail_bound=tails$bound)
}

cdf_survival <- function(c,a,n=64,L=-48,U=6) {
 z<-analytic(c,a,'survival',0);if(!is.null(z))return(z)
 cuts<-sort(unique(c(L,U,features[features>L&features<U],
  c+a*c(-1000,-100,-10,-1,0,1,10,100,1000))))
 cuts<-cuts[cuts>=L&cuts<=U];q<-gl_rule(n);half<-diff(cuts)/2;mid<-cuts[-length(cuts)]+half
 z<-sweep(outer(q$nodes,half,'*'),2,mid,'+')
 lw<-sweep(matrix(log(q$weights),nrow=n,ncol=length(half)),2,log(half),'+')
 lv<-lsum(c(z-exp(z)+pcauchy(z,location=c,scale=a,log.p=TRUE)+lw))
 bound<-exp(L)+exp(-exp(U))
 finish(exp(lv),messages='OK_CDF_IDENTITY_FIXED_RULE',evaluations=length(z),
  logvalue=lv,tail_bound=bound)
}

