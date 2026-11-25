# Selectively extracted definitions; no historical runner is sourced.
# See docs/PROVENANCE.csv for source identity and extraction scope.
A03a_roster <- function() {
  data.frame(setting = c("null", "alternative"), family = "Same",
             trial_N = 300L, r = 2, external_sizes = "100,250,400",
             delta = c(0, -0.5), generating_conditional_HR = exp(c(0, -0.5)),
             planned_replications = 1000L)
}

A03a_config <- function() {
  list(roster = A03a_roster(), external_sizes = c(100L, 250L, 400L),
       trial_beta = c(.20, .30, .40),
       external_beta = rbind(c(.20, .30, .40), c(.15, .35, .45), c(.22, .28, .38)),
       lambda = .01, nu = 1.2, total_time = 60, dropout = .15,
       times_grid = 1:60, r = 2, trial_N = 300L,
       n_calib = 100L, calib_steps = 1000L, calib_burnin = 200L,
       calib_pred_draws = 50L, final_steps = 5000L, final_burnin = 1000L,
       n_thin_draws = 1000L, init_par = c(1, -2, 0, -4.5),
       production_authorized = FALSE, schema = "A03a-paired-v2")
}

A03a_patient_summaries <- function(draws) {
  d <- dim(draws)
  stopifnot(length(d) == 3L, all(d > 0))
  mat <- matrix(draws, nrow = d[1] * d[2], ncol = d[3])
  list(median = matrix(matrixStats::rowMedians(mat), d[1], d[2]))
}

A03a_calibrate <- function(pvalues) {
  # Preserve legacy type-7, na.rm=TRUE and cap, including its limitations.
  raw <- as.numeric(quantile(pvalues, .05, na.rm = TRUE, type = 7))
  threshold <- raw
  if (threshold > .05) threshold <- .05
  list(raw_threshold = raw, threshold = threshold,
       n_splits = length(pvalues), n_missing = sum(is.na(pvalues)),
       n_nonfinite = sum(!is.finite(pvalues)))
}

A03a_memory_usage <- function() {
  tab <- gc()
  # R builds may include a "limit (Mb)" column. The last column is the
  # memory-unit column paired with "max used", not a fixed column number.
  list(gc_table = tab, gc_max_used_MiB = sum(tab[, ncol(tab)]))
}

A03a_stage <- function(trial_data, source_coefs, train_ids, eval_ids,
                       config, engine, final = FALSE, tracker = NULL,
                       stage_id = if (final) "final" else "calibration") {
  stage_start <- proc.time()[["elapsed"]]
  methods <- "median"
  n_source <- length(source_coefs)
  expert_curves <- setNames(lapply(methods, function(m)
    matrix(NA_real_, length(config$times_grid), n_source)), methods)
  gammas <- numeric(n_source)
  fit_diagnostics <- vector("list", n_source)
  for (s in seq_len(n_source)) {
    fit_start <- proc.time()[["elapsed"]]
    record_id <- paste(stage_id, s, sep = "/source_")
    if (!is.null(tracker)) {
      tracker$current <- record_id
      tracker$fit_records[[record_id]] <- list(status = "started")
    }
    f_train <- as.numeric(trial_data$X_mat[train_ids, , drop = FALSE] %*%
                            source_coefs[[s]])
    fit <- engine$mcmc_survival_ph(
      y = trial_data$time[train_ids], cens = trial_data$status[train_ids],
      f = f_train, init_par = config$init_par,
      steps = if (final) config$final_steps else config$calib_steps,
      burnin = if (final) config$final_burnin else config$calib_burnin)
    fit_elapsed <- proc.time()[["elapsed"]] - fit_start
    fit_finite <- all(is.finite(fit$chain))
    if (!is.null(tracker)) tracker$fit_records[[record_id]] <- list(
      status = if (fit_finite) "returned_finite_chain" else "returned_nonfinite_chain",
      elapsed_seconds = fit_elapsed, n_chain_rows = nrow(fit$chain))
    gammas[s] <- median(exp(fit$chain[, 2]))
    n_draws <- nrow(fit$chain)
    keep <- if (final && n_draws <= config$n_thin_draws) seq_len(n_draws) else
      round(seq(1, n_draws, length.out = if (final)
        config$n_thin_draws else config$calib_pred_draws))
    f_eval <- as.numeric(trial_data$X_mat[eval_ids, , drop = FALSE] %*%
                           source_coefs[[s]])
    f_eps <- 1e-8
    small <- abs(f_eval) < f_eps
    f_eval[small] <- sign(f_eval[small]) * f_eps
    prediction_start <- proc.time()[["elapsed"]]
    draws <- engine$posterior_predict_srv_prob_ph_optimized(
      fit$chain[keep, , drop = FALSE], f_eval, config$times_grid, n_rep = 1)
    draw_checks <- list(finite = all(is.finite(draws)),
      in_unit_interval = all(is.finite(draws) & draws >= 0 & draws <= 1),
      nonincreasing = all(draws[, -1, , drop = FALSE] <=
                           draws[, -length(config$times_grid), , drop = FALSE] + 1e-12))
    summaries <- A03a_patient_summaries(draws)
    for (method in methods)
      expert_curves[[method]][, s] <- colMeans(summaries[[method]])
    fit_diagnostics[[s]] <- list(n_chain_rows = n_draws, prediction_rows = keep,
      final_accept_rate = fit$final_accept_rate,
      fitting_status = if (fit_finite) "returned_finite_chain" else "returned_nonfinite_chain",
      fit_elapsed_seconds = fit_elapsed,
      prediction_elapsed_seconds = proc.time()[["elapsed"]] - prediction_start,
      conditional_survival_checks = draw_checks)
  }
  raw_weights <- 1 / (gammas + 1e-6)
  weights <- raw_weights / sum(raw_weights)
  results <- setNames(lapply(methods, function(method) {
    curve <- as.vector(expert_curves[[method]] %*% weights)
    stats <- engine$compare_predicted_to_observed(
      curve, config$times_grid, trial_data$time[eval_ids],
      trial_data$status[eval_ids])
    stats$log_HR <- log(stats$HR)
    stats$expected_by_patient <- approx(config$times_grid, -log(curve),
                                        xout = trial_data$time[eval_ids], rule = 2)$y
    stats$cohort_survival <- curve
    stats$source_cohort_survival <- expert_curves[[method]]
    stats$indicators <- list(zero_events = stats$O == 0,
      nonfinite_result = any(!is.finite(unlist(stats[c("O", "E", "HR", "Z", "p_value")]))),
      nonfinite_log_HR = !is.finite(stats$log_HR),
      positive_finite_E = is.finite(stats$E) && stats$E > 0,
      survival_in_unit_interval = all(is.finite(curve) & curve >= 0 & curve <= 1),
      survival_nonincreasing = all(diff(curve) <= 1e-12))
    stats
  }), methods)
  ans <- list(stage_id = stage_id, train_ids = train_ids, eval_ids = eval_ids,
       obs_time = trial_data$time[eval_ids], obs_status = trial_data$status[eval_ids],
       source_scales = gammas, source_weights = weights,
       fit_diagnostics = fit_diagnostics, methods = results,
       status = "completed", elapsed_seconds = proc.time()[["elapsed"]] - stage_start)
  if (!is.null(tracker)) tracker$stages[[stage_id]] <- ans
  ans
}

A03a_analyze_pair <- function(trial_data, source_coefs, config, engine, tracker = NULL) {
  control_ids <- trial_data$data$id[trial_data$data$trt == 0]
  treated_ids <- trial_data$data$id[trial_data$data$trt == 1]
  stopifnot(identical(trial_data$data$id, seq_len(nrow(trial_data$data))),
            length(control_ids) > 0L, length(treated_ids) > 0L)
  split_results <- vector("list", config$n_calib)
  for (b in seq_len(config$n_calib)) {
    train_indices <- sample(seq_along(control_ids), floor(length(control_ids) /
                                                          (config$r + 1)))
    test_indices <- setdiff(seq_along(control_ids), train_indices)
    split_results[[b]] <- A03a_stage(trial_data, source_coefs,
      control_ids[train_indices], control_ids[test_indices], config, engine,
      tracker = tracker, stage_id = sprintf("split_%03d", b))
    if (!is.null(tracker) && b %% 10L == 0L)
      message("Completed calibration split ", b, "/", config$n_calib)
  }
  final <- A03a_stage(trial_data, source_coefs, control_ids, treated_ids,
                     config, engine, final = TRUE, tracker = tracker)
  calibration <- setNames(lapply("median", function(method) {
    pvalues <- vapply(split_results, function(x) x$methods[[method]]$p_value,
                      numeric(1))
    ans <- A03a_calibrate(pvalues)
    ans$pvalues <- pvalues
    ans$reject <- final$methods[[method]]$p_value < ans$threshold
    ans
  }), "median")
  list(schema = config$schema, final = final, calibration = calibration,
       splits = split_results, source_coefs = source_coefs)
}

