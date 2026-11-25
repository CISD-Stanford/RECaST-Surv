# Later author-approved ALS reporting; kept separate from the original lower tail.
als_release_rules <- function(p,splits) {
  stopifnot(length(p)==1L,length(splits)==100L,all(is.finite(c(p,splits))),all(c(p,splits)>=0&c(p,splits)<=1))
  original_cutoff<-min(as.numeric(quantile(splits,.05,type=7,na.rm=FALSE)),.05)
  two_final<-2*min(p,1-p);two_splits<-2*pmin(splits,1-splits)
  cutoff_C<-min(as.numeric(quantile(two_splits,.05,type=7,na.rm=FALSE)),.05)
  list(original_lower_tail=list(p=p,cutoff=original_cutoff,reject=p<original_cutoff,uncalibrated_reject=p<.05),
    hybrid_Option_C=list(p=two_final,cutoff=cutoff_C,reject=two_final<cutoff_C,uncalibrated_reject=two_final<.05))
}
als_release_two_sided_cox <- function(result) {
  result$original_lower_tail_p<-result$p_value
  result$p_value<-2*pnorm(-abs(result$beta/result$SE));result$reject<-result$p_value<.05
  result$reporting_rule<-'Author-approved two-sided Wald amendment, 2026-11-20'
  result
}
