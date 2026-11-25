# Selectively extracted definitions; no historical runner is sourced.
# See docs/PROVENANCE.csv for source identity and extraction scope.
A07_ci <- function(O,E,sO,sE,B) {
  stopifnot(B %in% A07_B)
  if(B==100)return(construct_ci(O,E,sO,sE))
  # Only the required bank size is generalized. A05 rules and inversion are reused.
  reasons <- character();q<-c(NA_real_,NA_real_);T<-rep(NA_real_,length(sO))
  if(length(sO)!=B||length(sE)!=B)reasons<-c(reasons,paste0("requires_exactly_",B,"_splits"))
  if(!is.numeric(sO)||any(!is.finite(sO))||any(sO<0)||any(sO!=floor(sO)))reasons<-c(reasons,"invalid_split_O")
  if(!is.numeric(sE)||any(!is.finite(sE))||any(sE<=0))reasons<-c(reasons,"nonpositive_or_nonfinite_split_E")
  if(!length(reasons)) {
    T<-(sO-sE)/sqrt(sE)
    if(!all(is.finite(T)))reasons<-c(reasons,"nonfinite_split_pivot")
    else q<-unname(quantile(T,c(.025,.975),type=7,na.rm=FALSE))
  }
  out<-invert_cutoffs(O,E,q[1],q[2])
  if(length(reasons)) {
    out$status<-"unavailable"
    out$reason<-paste(unique(c(reasons,strsplit(out$reason,";",fixed=TRUE)[[1]])),collapse=";")
  }
  c(out,list(qL=q[1],qU=q[2],q_sum=sum(q),finite_split_pivots=sum(is.finite(T)),T=list(T)))
}

A07_calibrate <- function(O,E,sO,sE,B,full_p) {
  ci <- A07_ci(O,E,sO,sE,B)
  T <- ci$T[[1]]
  valid <- length(T)==B && all(is.finite(T)) && length(sO)==B && length(sE)==B &&
    all(is.finite(sO)) && all(sO>=0 & sO==floor(sO)) && all(is.finite(sE)&sE>0)
  p <- if(valid)pnorm(T) else rep(NA_real_,length(T))
  raw <- if(valid)as.numeric(quantile(p,.05,type=7,na.rm=FALSE)) else NA_real_
  cap <- if(valid)min(raw,.05) else NA_real_
  list(B=B,raw_threshold=raw,capped_threshold=cap,
    threshold_status=if(valid)"available" else "unavailable",
    reject=if(valid && valid_count(O) && valid_expected(E) && length(full_p)==1 && is.finite(full_p))full_p<cap else NA,
    ci=ci,n_planned=B,n_observed=length(sO),n_finite=sum(is.finite(T)),
    n_zero=if(all(is.finite(sO)))sum(sO==0) else NA_integer_)
}

