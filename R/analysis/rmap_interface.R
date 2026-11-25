# Selectively extracted definitions; no historical runner is sourced.
# See docs/PROVENANCE.csv for source identity and extraction scope.
convert_ind_data_to_poi_data <- function(inddata){
  r <- sum(inddata$status == 1, na.rm = T)
  e <- sum(inddata$time, na.rm = T)
  return(list(r = r, e = e))
}

rMAP_surv_fit = function(mydata, prior.std = 2, w = 0.5 , parameters = "theta.new", niter = 5000) {
  
  n.hist = length(mydata$r)-2
  
  dataTemp = list(
    "r" = mydata$r,
    "E" = mydata$e,
    "n.hist" = n.hist,
    "r.trt" = mydata$r[length(mydata$r)],
    "E.trt" = mydata$e[length(mydata$e)],
    'std_heterogeneity' = prior.std,
    "weights" = c(w,1-w),
    "scaling.factor" = 10
  )
  
  model <-
    jags.model(
      file = "rMAP_Surv_Pois_Fit_ctl.txt",
      data = dataTemp,
      n.chains = 1,
      n.adapt = 0.2 * niter,
      quiet = TRUE
    )
  update(model, n.iter = 0.2 * niter, progress.bar = 'none') # burn in
  MAP_model <-
    coda.samples(model,
                 variable.names = parameters,
                 thin = 5,
                 n.iter = niter, progress.bar = 'none')
  
  return(MAP_model)
}

params <- c("beta1", "log_lambda.pred", "theta")

