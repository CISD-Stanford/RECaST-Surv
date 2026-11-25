# Selectively extracted definitions; no historical runner is sourced.
# See docs/PROVENANCE.csv for source identity and extraction scope.
AD_generate <- function(engine, cfg) {
  bias_config <- function(s) switch(s, "ALL"=NULL, "[0,Inf)"=list(var=1, range=c(0,Inf)),
                                    "(-Inf,0]"=list(var=1, range=c(-Inf,0)), stop("Unknown X1 support"))
  generate <- function(g) {
    is_trial <- g==1L
    j <- g-1L
    rho <- cfg$AD$rho[g]
    args <- list(n=if(is_trial) cfg$trial_N else cfg$external_sizes[j],
      source=if(is_trial) "trial" else "external", ext_model_type="same",
      lambda=cfg$lambda, nu=cfg$nu, r=cfg$r, delta=cfg$delta,
      betas=if(is_trial) cfg$trial_beta else cfg$external_beta[j,],
      total_time=cfg$total_time, random_cens_rate=-log(1-cfg$dropout)/cfg$total_time,
      cor_list=if(rho==0) NA else list(c(1,2,rho)))
    if (startsWith(cfg$AD$case, "D")) args["bias_config"] <- list(bias_config(cfg$AD$support[g]))
    do.call(engine$simulate_survival_data, args)
  }
  trial <- generate(1L)
  list(trial=trial, external=lapply(2:4, generate))
}

AD_source_coefficients <- function(models, cfg) {
  raw <- lapply(models, stats::coef)
  stopifnot(length(raw)==3L, all(lengths(raw)==3L), all(is.finite(unlist(raw))))
  lapply(seq_along(raw), function(j) {
    b <- raw[[j]]
    settings <- cfg$AD$corruption[[as.character(j)]]
    if (is.null(settings)) return(b)
    if (settings$bias!=0) b <- b+settings$bias
    if (settings$noise_sd>0) b <- b+stats::rnorm(length(b), mean=0, sd=settings$noise_sd)
    b
  })
}

