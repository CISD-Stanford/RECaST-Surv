# Selectively extracted definitions; no historical runner is sourced.
# See docs/PROVENANCE.csv for source identity and extraction scope.
decision <- function(p) p<.05

wald_fields <- function(coefficients,V) {
  stopifnot('treatment' %in% names(coefficients),'treatment' %in% rownames(V),
            identical(rownames(V),colnames(V)))
  beta <- unname(coefficients['treatment']); variance <- unname(V['treatment','treatment']); se <- sqrt(variance)
  stopifnot(is.finite(beta),is.finite(se),se>0)
  z <- beta/se; W <- z*z; p <- 2*stats::pnorm(-abs(z)); p2 <- stats::pchisq(W,1,lower.tail=FALSE)
  list(beta=beta,model_se=se,treatment_variance=variance,HR=exp(beta),z=z,W=W,
       p=p,p_chisq=p2,reject=decision(p),lower_log=beta-stats::qnorm(.975)*se,
       upper_log=beta+stats::qnorm(.975)*se,CI_lower=exp(beta-stats::qnorm(.975)*se),
       CI_upper=exp(beta+stats::qnorm(.975)*se))
}

fit_once <- function(trial,method) {
  dat <- trial$data
  formula <- if(method=='Cox-U') survival::Surv(time,event) ~ treatment else
    survival::Surv(time,event) ~ treatment+X1+X2+X3
  warnings <- character(); error <- NULL; initial_rng <- .Random.seed
  fit <- tryCatch(withCallingHandlers(survival::coxph(formula,data=dat,ties='efron',
            robust=FALSE,na.action=stats::na.fail,singular.ok=TRUE,
            control=survival::coxph.control(eps=1e-9,iter.max=20,timefix=TRUE),
            model=FALSE,x=FALSE,y=FALSE),warning=function(w) {
              warnings <<- c(warnings,conditionMessage(w)); invokeRestart('muffleWarning')
            }),error=function(e) { error <<- conditionMessage(e); NULL })
  stopifnot(identical(initial_rng,.Random.seed))
  out <- list(schema='COX_ALL_FIT_v001',trial_id=trial$trial_id,bank_id=trial$bank_id,
       replication_id=trial$replication_id,method=method,data_hash=trial$data_hash,
       source_identity=trial$source_identity,formula=paste(deparse(formula),collapse=' '),
       ties='efron',covariance_type='ordinary_model_based',robust=FALSE,attempted=TRUE,
       warnings=warnings,error=error,status='NON_ESTIMABLE',inference_available=FALSE,
       coefficients=NULL,named_vcov=NULL,summary_coefficients=NULL,fit_metadata=NULL)
  reasons <- character()
  if(is.null(fit)) reasons <- c(reasons,paste0('FIT_ERROR: ',error)) else {
    co <- stats::coef(fit); V <- stats::vcov(fit); sm <- summary(fit)
    out$coefficients <- co; out$named_vcov <- V; out$summary_coefficients <- sm$coefficients
    out$fit_metadata <- list(n=fit$n,nevent=fit$nevent,iter=fit$iter,loglik=fit$loglik,
        score=fit$score,method=fit$method,has_naive_var=!is.null(fit$naive.var),
        summary_used_robust=isTRUE(sm$used.robust),call=paste(deparse(fit$call),collapse=' '))
    expected <- if(method=='Cox-U') 'treatment' else c('treatment','X1','X2','X3')
    if(!identical(names(co),expected)) stop('STRUCTURAL: coefficient names mismatch')
    if(!all(is.finite(co))) reasons <- c(reasons,'NONFINITE_COEFFICIENT')
    if(!all(is.finite(V)) || any(diag(V)<=0) || !isTRUE(all.equal(V,t(V),tolerance=1e-12))) reasons <- c(reasons,'INVALID_VARIANCE')
    if(!is.null(fit$naive.var) || isTRUE(sm$used.robust)) stop('STRUCTURAL: robust covariance found')
    if(any(grepl('converg|infinite|infinity|ran out|singular|iteration',warnings,ignore.case=TRUE))) reasons <- c(reasons,'CONVERGENCE_OR_INFINITY_WARNING')
    if(!length(reasons)) {
      ww <- tryCatch(wald_fields(co,V),error=function(e) NULL)
      if(is.null(ww)) reasons <- c(reasons,'INVALID_TREATMENT_SE') else {
        stopifnot(identical(ww$model_se,unname(sm$coefficients['treatment','se(coef)'])))
        stopifnot(abs(ww$p-ww$p_chisq)<5e-14)
        out <- c(out,ww); out$inference_available <- TRUE
        out$status <- if(length(warnings)) 'ESTIMABLE_WITH_WARNING' else 'ESTIMABLE'
      }
    }
  }
  out$failure_reasons <- reasons
  out
}

