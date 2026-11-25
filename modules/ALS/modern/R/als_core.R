# In-memory library. All data-bearing orchestration requires an explicitly
# synthetic context. The context is an execution boundary, not data authentication.
als_assert <- function(ok, message) {
  if (length(ok) != 1L || is.na(ok) || !ok) stop(message, call. = FALSE)
}

als_synthetic_context <- function(fixture_id) {
  als_assert(is.character(fixture_id) && length(fixture_id) == 1L && nzchar(fixture_id),
             "A named synthetic fixture is required")
  structure(list(origin = "synthetic", fixture_id = fixture_id), class = "als_context")
}

als_require_synthetic <- function(context) {
  als_assert(inherits(context, "als_context") && identical(context$origin, "synthetic"),
    "STOP: v001 permits synthetic fixtures only; real-data execution is not authorized")
}

als_outcomes <- function(time, status, horizon = 120) {
  als_assert(is.numeric(time) && length(time) > 0L && length(time) == length(status) &&
    all(is.finite(time)) && all(time > 0) && all(!is.na(status) & status %in% c(0, 1)),
    "Invalid observed time/status; no eligibility filtering or imputation is permitted")
  als_assert(identical(as.numeric(horizon), 120), "Only the approved 120-month horizon is supported")
  data.frame(time = pmin(time, horizon), status = as.integer(status * (time <= horizon)),
             administratively_truncated = time > horizon,
             outcome_origin = "observed_registry", censoring_origin = "observed_plus_administrative")
}

als_standardize <- function(df, context, config = als_config()) {
  als_require_synthetic(context)
  als_assert(is.data.frame(df) && !anyDuplicated(names(df)) && "SEX" %in% names(df),
             "Raw SEX and unambiguous columns are required")
  als_assert(all(!is.na(df$SEX) & df$SEX %in% c("M", "F")),
             "Unexpected/missing raw SEX category; only exact M and F are allowed")
  sex <- as.integer(df$SEX == "M")
  if ("SEX_NUM" %in% names(df))
    als_assert(identical(as.numeric(df$SEX_NUM), as.numeric(sex)), "SEX_NUM conflicts with raw SEX")
  df$SEX_NUM <- sex
  covs <- config$covariates
  als_assert(all(covs %in% names(df)), "Missing approved covariate")
  als_assert(all(vapply(df[covs], is.numeric, logical(1))), "Covariates must be numeric")
  raw_X <- as.matrix(df[covs])
  als_assert(all(is.finite(raw_X)), "Incomplete/nonfinite covariates; STOP rather than drop rows")
  centers <- colMeans(raw_X)
  scales <- apply(raw_X, 2, sd)
  als_assert(all(is.finite(scales) & scales > 0), "Degenerate covariate; no SD substitution")
  X <- sweep(sweep(raw_X, 2, centers, "-"), 2, scales, "/")
  list(X = X, raw_X = raw_X, centers = centers, sample_sd = scales,
       covariates = covs, n = nrow(X), context = context)
}

als_prepare_fixture <- function(df, patient_ids, context, config = als_config()) {
  als_require_synthetic(context)
  als_assert(nrow(df) == config$n_registry && length(patient_ids) == nrow(df),
             "Exactly 2790 eligible records and identity keys are required; no surplus")
  als_assert(!anyNA(patient_ids) && !anyDuplicated(patient_ids) &&
               all(nzchar(as.character(patient_ids))), "Identity/uniqueness assumption failed")
  p <- als_standardize(df, context, config)
  als_assert(all(c("STIME", "STATUS") %in% names(df)), "Missing STIME/STATUS")
  p$outcomes <- als_outcomes(df$STIME, df$STATUS, config$horizon)
  p$row_ids <- seq_len(nrow(df))
  p$patient_ids <- patient_ids
  p
}

# Cox cumulative baseline is a right-continuous step function. The generalized
# inverse is inf{t: H0(t) >= target}, not linear interpolation of hazard jumps.
als_validate_baseline <- function(baseline, horizon = 120) {
  als_assert(is.data.frame(baseline) && all(c("time", "hazard") %in% names(baseline)) &&
    nrow(baseline) > 0L, "Missing baseline table")
  als_assert(all(is.finite(baseline$time)) && all(is.finite(baseline$hazard)) &&
    all(baseline$time >= 0 & baseline$time <= horizon) &&
    all(diff(baseline$time) > 0) && all(baseline$hazard >= 0) &&
    all(diff(baseline$hazard) >= 0), "Invalid baseline support or cumulative hazards")
  invisible(TRUE)
}

als_invert_baseline <- function(target_H, baseline, horizon = 120) {
  als_validate_baseline(baseline, horizon)
  als_assert(is.numeric(target_H) && all(!is.na(target_H) & target_H >= 0), "Invalid target hazard")
  H120 <- tail(baseline$hazard, 1L) # flat after last supported jump, through 120
  vapply(target_H, function(h) {
    if (h == 0) return(0)
    if (h > H120) return(Inf)
    baseline$time[which(baseline$hazard >= h)[1L]]
  }, numeric(1))
}

als_new_generator_state <- function() {
  e <- new.env(parent = emptyenv())
  e$fit_count <- 0L
  e$attempted <- FALSE
  e
}

als_capture_fit <- function(expr) {
  warnings <- character()
  fit <- withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  # No fit-warning acceptance/retry policy was approved. Stop this branch.
  als_assert(length(warnings) == 0L, paste("STOP: Cox warning; no approved retry policy:",
                                         paste(warnings, collapse = "; ")))
  fit
}

als_fit_generator_once <- function(state, prepared, context, config = als_config()) {
  als_require_synthetic(context)
  als_assert(prepared$n == config$n_registry && identical(colnames(prepared$X), config$covariates),
             "Generator requires all 2790 rows and exactly seven covariates")
  if (state$fit_count > 0L) {
    als_assert(identical(state$prepared, prepared), "Fixed generator cannot be refitted to changed data")
    return(state$cache)
  }
  als_assert(!state$attempted, "STOP: global generator attempt already failed; no retry policy is approved")
  state$attempted <- TRUE
  d <- data.frame(time = prepared$outcomes$time, status = prepared$outcomes$status, prepared$X)
  formula <- stats::as.formula(paste("survival::Surv(time,status) ~", paste(config$covariates, collapse = "+")))
  model <- als_capture_fit(survival::coxph(formula, data = d, ties = "efron",
    na.action = stats::na.fail, singular.ok = FALSE, x = TRUE, y = TRUE, model = TRUE))
  beta <- stats::coef(model)[config$covariates]
  als_assert(model$n == config$n_registry && all(is.finite(beta)), "Incomplete/invalid global generator")
  baseline <- survival::basehaz(model, centered = FALSE)
  als_validate_baseline(baseline, config$horizon)
  cache <- list(beta = beta, baseline = baseline, H120 = tail(baseline$hazard, 1L),
    model = model, centered = FALSE, centers = prepared$centers, sample_sd = prepared$sample_sd,
    n = model$n, horizon = config$horizon, role = "fixed_empirical_working_DGM")
  state$prepared <- prepared
  state$cache <- cache
  state$fit_count <- 1L
  cache
}

als_generate_treated <- function(X, generator, delta, U, C_random, context, config = als_config()) {
  als_require_synthetic(context)
  als_assert(is.matrix(X) && identical(colnames(X), config$covariates) &&
    identical(names(generator$beta), config$covariates) && all(is.finite(X)) &&
    all(is.finite(generator$beta)) && identical(generator$centered, FALSE),
    "Generator must use zero-reference baseline and exactly aligned standardized X beta")
  als_assert(length(delta) == 1L && delta %in% c(-0.5, 0), "Unapproved Delta")
  n <- nrow(X)
  als_assert(length(U) == n && all(is.finite(U) & U > 0 & U < 1) &&
    length(C_random) == n && all(!is.na(C_random) & C_random >= 0), "Invalid event/censoring draws")
  lp <- drop(X %*% generator$beta) # no second subtraction of coxph$means
  als_assert(all(is.finite(lp + delta)), "Nonfinite generating linear predictor")
  risk <- exp(lp + delta)
  als_assert(all(is.finite(risk) & risk > 0), "Numerical generating-risk overflow/underflow; no repair")
  target_H <- -log(U) / risk
  als_assert(all(!is.na(target_H) & target_H > 0), "Numerical target underflow; no time-zero event repair")
  T_event <- als_invert_baseline(target_H, generator$baseline, config$horizon)
  C <- pmin(C_random, config$horizon)
  event <- is.finite(T_event) & T_event <= C
  data.frame(time = pmin(T_event, C), status = as.integer(event), T_event = T_event,
    C_random = C_random, C = C, target_H = target_H,
    administratively_truncated = !event & C_random >= config$horizon,
    outcome_origin = "generated_treated", censoring_origin = "independent_exponential_plus_administrative")
}

als_trial_weights <- function(X, config = als_config()) {
  als_assert(is.matrix(X) && all(c("AGE", "DISDUR") %in% colnames(X)) && all(is.finite(X)),
             "Full-registry standardized AGE/DISDUR are required")
  ell <- config$age_log_weight * X[, "AGE"] + config$disdur_log_weight * X[, "DISDUR"]
  w <- exp(ell - max(ell))
  cap <- config$weight_cap_median_multiple * median(w)
  capped <- pmin(w, cap)
  als_assert(is.finite(cap) && cap > 0 && all(is.finite(capped)) && sum(capped) > 0 &&
    sum(capped > 0) >= config$n_trial, "Weight underflow/availability failure; no weight repair")
  list(ell = ell, uncapped = w, cap = cap, capped = capped, probability = capped / sum(capped))
}

als_validate_assignments <- function(a, n = 2790L, config = als_config()) {
  als_assert(n == config$n_registry && n == config$n_trial + sum(config$n_external),
             "Allocation has no surplus; eligibility count must be exactly 2790")
  als_assert(length(a$trial) == config$n_trial && !anyDuplicated(a$trial) &&
    length(a$trt) == config$n_trial && all(a$trt %in% c(0L, 1L)) &&
    sum(a$trt == 1L) == config$n_treated && sum(a$trt == 0L) == config$n_control,
    "Invalid trial size or exact treatment allocation")
  als_assert(identical(names(a$external), names(config$n_external)) &&
    identical(as.integer(lengths(a$external)), as.integer(config$n_external)), "External sizes/names differ")
  all_ids <- c(a$trial, unlist(a$external, use.names = FALSE))
  als_assert(is.numeric(all_ids) && all(is.finite(all_ids)) && !anyDuplicated(all_ids) &&
    setequal(all_ids, seq_len(n)), "Cohorts must be disjoint and exhaust all available rows")
  invisible(TRUE)
}

als_sample_assignments <- function(prepared, design, context, config = als_config()) {
  als_require_synthetic(context)
  als_assert(design %in% c("random", "biased") && length(design) == 1L, "Unknown design")
  als_assert(prepared$n == config$n_registry && !anyDuplicated(prepared$patient_ids) &&
    length(prepared$patient_ids) == config$n_registry && !anyNA(prepared$patient_ids),
    "Eligibility/unique-patient availability failed")
  weights <- if (design == "biased") als_trial_weights(prepared$X, config) else NULL
  trial <- sample.int(prepared$n, config$n_trial, replace = FALSE,
                      prob = if (is.null(weights)) NULL else weights$capped)
  trt <- integer(config$n_trial)
  trt[sample.int(config$n_trial, config$n_treated, replace = FALSE)] <- 1L
  remaining <- setdiff(seq_len(prepared$n), trial)
  external <- lapply(config$n_external, function(size) {
    # sample.int avoids R's sample(length-one-numeric) special case.
    chosen <- remaining[sample.int(length(remaining), size, replace = FALSE)]
    remaining <<- setdiff(remaining, chosen)
    chosen
  })
  a <- list(design = design, trial = trial, trt = trt, external = external,
            weight_details = weights, fixture_id = context$fixture_id)
  als_validate_assignments(a, prepared$n, config)
  a
}

als_paired_scenarios <- function(prepared, assignments, generator, context, draws = NULL,
                                 config = als_config()) {
  als_require_synthetic(context)
  als_validate_assignments(assignments, prepared$n, config)
  treated <- which(assignments$trt == 1L)
  if (is.null(draws)) {
    # Successive independent RNG calls; no resetting seeds or shared U/C draws.
    draws <- lapply(c("ALTERNATIVE", "DIAGNOSTIC"), function(s)
      list(U = stats::runif(config$n_treated), C_random = stats::rexp(config$n_treated, config$censor_rate)))
    names(draws) <- c("ALTERNATIVE", "DIAGNOSTIC")
  }
  als_assert(identical(sort(names(draws)), sort(c("ALTERNATIVE", "DIAGNOSTIC"))),
             "Supply separate explicitly named alternative and diagnostic draws")
  scenarios <- lapply(config$scenarios, function(s) {
    trial <- prepared$outcomes[assignments$trial, , drop = FALSE]
    trial$row_id <- assignments$trial
    trial$trt <- assignments$trt
    trial$T_event <- trial$C_random <- trial$C <- trial$target_H <- NA_real_
    if (s != "NULL") {
      g <- als_generate_treated(prepared$X[assignments$trial[treated], , drop = FALSE], generator,
        config$delta[[s]], draws[[s]]$U, draws[[s]]$C_random, context, config)
      trial[treated, names(g)] <- g
    }
    trial$scenario <- s
    external <- lapply(assignments$external, function(ids) {
      d <- prepared$outcomes[ids, , drop = FALSE]; d$row_id <- ids; d$scenario <- s; d
    })
    list(scenario = s, assignments = assignments, trial = trial, external = external,
         interpretation = als_metric_labels(s))
  })
  stats::setNames(scenarios, config$scenarios)
}

als_source_weights <- function(scales, epsilon = 1e-6) {
  als_assert(length(scales) == 3L && all(is.finite(scales) & scales > 0), "Invalid source scales")
  raw <- 1 / (scales + epsilon)
  raw / sum(raw)
}

als_aggregate_predictions <- function(draw_arrays, scales, config = als_config()) {
  als_assert(length(draw_arrays) == 3L, "Exactly three source prediction arrays are required")
  expected_dim <- dim(draw_arrays[[1L]])[1:2]
  curves <- vapply(draw_arrays, function(a) {
    als_assert(length(dim(a)) == 3L && all(dim(a) > 0) &&
      identical(dim(a)[1:2], expected_dim) && dim(a)[2L] == length(config$time_grid) &&
      all(is.finite(a) & a >= 0 & a <= 1), "Invalid [patient,time,draw] predictions")
    d <- dim(a)
    dim(a) <- c(d[1L] * d[2L], d[3L])
    als_assert(requireNamespace("matrixStats", quietly = TRUE), "matrixStats is required; no package installation")
    patient_medians <- matrix(matrixStats::rowMedians(a), nrow = d[1L], ncol = d[2L])
    als_assert(all(apply(patient_medians, 1L, function(x) all(diff(x) <= 1e-12))),
               "Survival curves must be nonincreasing")
    colMeans(patient_medians)
  }, numeric(length(config$time_grid)))
  weights <- als_source_weights(scales, config$source_weight_epsilon)
  combined <- as.vector(curves %*% weights)
  list(source_curves = curves, scales = scales, weights = weights, combined_survival = combined)
}

als_oe <- function(combined_survival, obs_time, obs_status, config = als_config()) {
  als_assert(length(combined_survival) == length(config$time_grid) &&
    all(is.finite(combined_survival) & combined_survival > 0 & combined_survival <= 1) &&
    all(diff(combined_survival) <= 1e-12), "Invalid survival curve; no clipping/pseudocount")
  als_assert(length(obs_time) > 0L && length(obs_time) == length(obs_status) &&
    all(is.finite(obs_time) & obs_time >= 0 & obs_time <= config$horizon) &&
    all(!is.na(obs_status) & obs_status %in% c(0L, 1L)), "Invalid evaluation outcomes")
  H <- -log(combined_survival)
  E_i <- stats::approx(config$time_grid, H, xout = obs_time, rule = 2)$y
  E <- sum(E_i); O <- sum(obs_status)
  HR <- O / E; Z <- (O - E) / sqrt(E); p <- stats::pnorm(Z)
  valid <- is.finite(HR) && HR > 0 && is.finite(E) && E > 0 && is.finite(p)
  list(O = O, E = E, E_individual = E_i, HR = HR, Z = Z, p_value = p,
       valid_HR = valid, status = if (valid) "ok" else "invalid_HR_or_expected_events")
}

als_calibration_threshold <- function(p_values, config = als_config()) {
  als_assert(length(p_values) == config$B && all(is.finite(p_values) & p_values >= 0 & p_values <= 1),
             "Exactly B=100 valid split p-values required; no NA removal")
  min(as.numeric(stats::quantile(p_values, probs = config$calibration_probability,
                                type = 7, na.rm = FALSE)), config$threshold_cap)
}

als_reject <- function(p_value, threshold) {
  als_assert(length(p_value) == 1L && is.finite(p_value) && p_value >= 0 && p_value <= 1 &&
    length(threshold) == 1L && is.finite(threshold) && threshold >= 0 && threshold <= 0.05,
    "Invalid final p-value/threshold")
  p_value < threshold
}

als_control_splits <- function(control_ids, context, config = als_config()) {
  als_require_synthetic(context)
  als_assert(length(control_ids) == 100L && !anyDuplicated(control_ids), "Exactly 100 concurrent controls required")
  lapply(seq_len(config$B), function(b) {
    train_pos <- sample.int(100L, config$calibration_train, replace = FALSE)
    list(split_id = b, train = control_ids[train_pos], held_out = control_ids[-train_pos])
  })
}

# Adapter contract: fit_target(y,cens,f,init_par,steps,burnin) returns list(chain),
# predict_draws(chain,f_treated,time_grid,n_rep) returns [patient,time,draw].
# The pinned audited engine is available through als_inherited_engine(), which
# requires a synthetic context. Tests use deterministic fit/prediction doubles.
als_recast_analysis <- function(trial, X, source_coefs, fit_target, predict_draws,
                                context, splits = NULL, config = als_config()) {
  als_require_synthetic(context)
  als_assert(nrow(trial) == 300L && nrow(X) == 300L && identical(colnames(X), config$covariates) &&
    sum(trial$trt == 0) == 100L && sum(trial$trt == 1) == 200L, "Invalid trial/X allocation")
  als_assert(length(source_coefs) == 3L && all(vapply(source_coefs, function(b)
    identical(names(b), config$covariates) && all(is.finite(b)), logical(1))), "Misaligned external coefficients")
  controls <- which(trial$trt == 0L); treated <- which(trial$trt == 1L)
  if (is.null(splits)) splits <- als_control_splits(controls, context, config)
  als_assert(length(splits) == 100L, "Exactly 100 control-only repeated splits required")
  stage <- function(train, evaluate, calibration) {
    arrays <- vector("list", 3L); scales <- numeric(3L)
    for (s in seq_len(3L)) {
      b <- source_coefs[[s]]
      chain <- fit_target(y = trial$time[train], cens = trial$status[train],
        f = drop(X[train, , drop = FALSE] %*% b), init_par = config$inherited_fit$init_par,
        steps = if (calibration) config$inherited_fit$calibration_steps else config$inherited_fit$final_steps,
        burnin = if (calibration) config$inherited_fit$calibration_burnin else config$inherited_fit$final_burnin)$chain
      als_assert(is.matrix(chain) && nrow(chain) > 0L && ncol(chain) == 4L && all(is.finite(chain)),
                 "Invalid target chain; no retry/drop policy")
      scales[s] <- median(exp(chain[, 2L])) # full retained chain, before draw subset
      n_keep <- if (calibration) config$inherited_fit$calibration_draws else config$inherited_fit$final_draws
      keep <- if (calibration || nrow(chain) > n_keep) round(seq(1, nrow(chain), length.out = n_keep)) else seq_len(nrow(chain))
      f <- drop(X[evaluate, , drop = FALSE] %*% b)
      f[abs(f) < config$score_epsilon] <- config$score_epsilon
      arrays[[s]] <- predict_draws(chain = chain[keep, , drop = FALSE], f_treated = f,
                                  time_grid = config$time_grid, n_rep = 1L)
      als_assert(identical(dim(arrays[[s]]), as.integer(c(length(evaluate), length(config$time_grid), length(keep)))),
                 "Predictor returned misaligned patient/time/draw dimensions")
    }
    a <- als_aggregate_predictions(arrays, scales, config)
    list(aggregation = a, oe = als_oe(a$combined_survival, trial$time[evaluate], trial$status[evaluate], config))
  }
  split_results <- lapply(splits, function(s) {
    als_assert(length(s$train) == 33L && length(s$held_out) == 67L &&
      !anyDuplicated(c(s$train, s$held_out)) && setequal(c(s$train, s$held_out), controls),
      "Split must partition concurrent controls only into 33/67")
    stage(s$train, s$held_out, TRUE)
  })
  p_values <- vapply(split_results, function(s) s$oe$p_value, numeric(1))
  threshold <- als_calibration_threshold(p_values, config)
  final <- stage(controls, treated, FALSE)
  list(method = "RECaST_Surv_OE_median", final = final, splits = splits,
       split_results = split_results, split_p_values = p_values, threshold = threshold,
       reject = als_reject(final$oe$p_value, threshold))
}

als_wald <- function(beta, se) {
  als_assert(length(beta) == 1L && is.finite(beta) && length(se) == 1L && is.finite(se) && se > 0,
             "Invalid treatment coefficient/SE")
  list(HR = exp(beta), beta = beta, SE = se, z = beta / se, p_value = stats::pnorm(beta / se))
}

als_cox_fit <- function(trial, X, adjusted, context, config = als_config()) {
  als_require_synthetic(context)
  als_assert(is.logical(adjusted) && length(adjusted) == 1L &&
    nrow(trial) == nrow(X) && identical(colnames(X), config$covariates) &&
    all(c("time", "status", "trt") %in% names(trial)), "Invalid comparator data")
  d <- data.frame(trial[c("time", "status", "trt")], X)
  terms <- if (adjusted) c("trt", config$covariates) else "trt"
  formula <- stats::as.formula(paste("survival::Surv(time,status) ~", paste(terms, collapse = "+")))
  model <- als_capture_fit(survival::coxph(formula, data = d, ties = "efron",
    na.action = stats::na.fail, singular.ok = FALSE, x = TRUE, y = TRUE, model = TRUE))
  result <- als_wald(stats::coef(model)[["trt"]], sqrt(stats::vcov(model)["trt", "trt"]))
  result$reject <- als_reject(result$p_value, config$cox_alpha)
  result$model <- model
  result$method <- if (adjusted) "Cox_trial_treatment_plus_seven_supplementary" else "Cox_trial_treatment_only_primary"
  result
}

als_fit_external <- function(outcomes, X, context, config = als_config()) {
  als_require_synthetic(context)
  als_assert(identical(colnames(X), config$covariates) && nrow(outcomes) == nrow(X),
             "External outcome/X mismatch")
  d <- data.frame(time = outcomes$time, status = outcomes$status, X)
  formula <- stats::as.formula(paste("survival::Surv(time,status) ~", paste(config$covariates, collapse = "+")))
  fit <- als_capture_fit(survival::coxph(formula, data = d, ties = "efron", na.action = stats::na.fail,
                                      singular.ok = FALSE, x = TRUE, y = TRUE, model = TRUE))
  als_assert(all(is.finite(stats::coef(fit))), "Invalid source coefficients")
  fit
}

# Bind only the five exact pinned definitions, never a historical script/runner.
# No fit runs at binding; only deterministic quadrature nodes are initialized.
als_inherited_engine <- function(implementation_root, context) {
  als_require_synthetic(context)
  for (package in c("statmod", "mvtnorm"))
    als_assert(requireNamespace(package, quietly = TRUE), paste("Missing required package", package))
  file <- file.path(implementation_root, "R", "inherited_engine_pinned.R")
  als_assert(file.exists(file) && identical(unname(tools::md5sum(file)), "0ddc2b6c4595a3b7c5f2fa341f421848"),
             "Pinned inherited-engine integrity check failed")
  e <- new.env(parent = asNamespace("stats"))
  sys.source(file, envir = e)
  list(fit_target = function(y, cens, f, init_par, steps, burnin) {
    als_require_synthetic(context)
    als_capture_fit(e$mcmc_survival_ph(y, cens, f, init_par, steps, burnin))
  }, predict_draws = function(chain, f_treated, time_grid, n_rep) {
    als_require_synthetic(context)
    als_assert(identical(as.integer(n_rep), 1L) && identical(as.integer(time_grid), 1:120),
               "Only approved n_rep=1 and grid1:120 are supported")
    e$posterior_predict_srv_prob_ph_optimized(chain, f_treated, time_grid, n_rep)
  }, environment = e, provenance = "Exact audited U5 Code/MCMC_PH.R selected definitions; inherited numerical limitations retained")
}
