# Load only release-local definitions. Loading performs no stochastic analysis.
release_load <- function(root) {
  e <- new.env(parent=globalenv())
  e$release_root <- normalizePath(root,mustWork=TRUE)
  e$coxph <- survival::coxph; e$Surv <- survival::Surv
  e$mvrnorm <- MASS::mvrnorm
  paths <- c('model/mcmc_ph.R','simulation/main.R','analysis/paired_prediction.R',
    'analysis/score_ci.R','analysis/cox.R','analysis/rmap.R','simulation/h_validate.R','simulation/h_dgm.R',
    'simulation/ad_adapter.R','simulation/dispatch.R','analysis/workflow.R',
    'sensitivity/calibration.R','sensitivity/workflow.R','reporting/summarize.R')
  for(p in paths)sys.source(file.path(e$release_root,'R',p),e)
  e$A07_B <- c(100L,200L,500L)
  e
}
