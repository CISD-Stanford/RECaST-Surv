# Reporting interfaces retain attempted, failed, invalid and missing records.
als_metric_labels <- function(scenario) {
  als_assert(length(scenario) == 1L && scenario %in% c("NULL", "ALTERNATIVE", "DIAGNOSTIC"),
             "Unknown scenario")
  switch(scenario,
    "NULL" = list(rejection = "null_type_I_error", squared_deviation = "null_MSE_relative_to_HR_1",
                  reference = 1, reference_role = "null_reference"),
    "ALTERNATIVE" = list(rejection = "alternative_rejection_probability_power_under_specified_generator",
      squared_deviation = "alternative_mean_squared_deviation_from_nominal_reference_NOT_true_parameter_MSE",
      reference = exp(-0.5), reference_role = "conditional_generator_multiplier_NOT_established_true_marginal_HR"),
    "DIAGNOSTIC" = list(rejection = "diagnostic_rejection_frequency_NOT_formal_type_I_error",
      squared_deviation = "not_defined_for_diagnostic", reference = NA_real_,
      reference_role = "hybrid_generated_observed_outcome_compatibility_NOT_genuine_null"))
}

als_result_record <- function(replication_id, design, scenario, method, HR = NA_real_,
                               p_value = NA_real_, threshold = NA_real_, status = "ok",
                               failure_reason = "", O = NA_real_, E = NA_real_,
                               split_valid_count = NA_integer_) {
  als_assert(scenario %in% c("NULL", "ALTERNATIVE", "DIAGNOSTIC") &&
    design %in% c("random", "biased") && length(replication_id) == 1L &&
    length(method) == 1L && nzchar(method), "Invalid result key")
  valid_HR <- length(HR) == 1L && is.finite(HR) && HR > 0
  valid_test <- length(p_value) == 1L && is.finite(p_value) && p_value >= 0 && p_value <= 1 &&
    length(threshold) == 1L && is.finite(threshold) && threshold >= 0 && threshold <= 0.05
  data.frame(replication_id = as.character(replication_id), design = design, scenario = scenario,
    method = method, HR = HR, p_value = p_value, threshold = threshold,
    reject = if (valid_test && status == "ok") p_value < threshold else NA,
    status = status, failure_reason = failure_reason, valid_HR = valid_HR,
    valid_test = valid_test, O = O, E = E, split_valid_count = split_valid_count)
}

# expected_ids is the approved roster, never inferred from successful outputs.
# Any missing/failed/nonpositive/nonfinite HR makes ALL full-roster HR summaries
# unavailable. No silent available-case estimate or pseudocount is emitted.
als_summarize_method <- function(records, expected_ids) {
  als_assert(nrow(records) > 0L && length(expected_ids) > 0L && !anyNA(expected_ids) &&
    !anyDuplicated(expected_ids) && !anyDuplicated(records$replication_id) &&
    all(records$replication_id %in% as.character(expected_ids)), "Invalid/duplicate reporting roster")
  for (v in c("design", "scenario", "method"))
    als_assert(length(unique(records[[v]])) == 1L, "Summarize one design/scenario/method at a time")
  labels <- als_metric_labels(records$scenario[1L])
  n_expected <- length(expected_ids)
  idx <- match(as.character(expected_ids), records$replication_id)
  r <- records[idx, , drop = FALSE]
  present <- !is.na(idx)
  fit_ok <- present & !is.na(r$status) & r$status == "ok"
  hr_ok <- fit_ok & is.finite(r$HR) & r$HR > 0
  test_ok <- fit_ok & is.finite(r$p_value) & r$p_value >= 0 & r$p_value <= 1 &
    is.finite(r$threshold) & r$threshold >= 0 & r$threshold <= 0.05
  all_hr <- all(hr_ok); all_test <- all(test_ok)
  h <- if (all_hr) r$HR else rep(NA_real_, n_expected)
  mse_terms <- if (all_hr && is.finite(labels$reference)) (h - labels$reference)^2 else rep(NA_real_, n_expected)
  reject_count <- sum((r$p_value < r$threshold)[test_ok])
  rate_available <- if (any(test_ok)) reject_count / sum(test_ok) else NA_real_
  metric <- function(x) if (all(is.finite(x))) mean(x) else NA_real_
  mcse <- function(x) if (length(x) > 1L && all(is.finite(x))) sd(x) / sqrt(length(x)) else NA_real_
  mean_log <- metric(log(h))
  geo <- exp(mean_log)
  list(design = records$design[1L], scenario = records$scenario[1L], method = records$method[1L],
    n_expected = n_expected, n_recorded = sum(present), n_missing = sum(!present),
    n_failed = sum(present & !fit_ok), n_HR_available = sum(hr_ok), n_test_available = sum(test_ok),
    n_nonfinite_HR = sum(present & !is.finite(r$HR)),
    n_nonpositive_HR = sum(present & is.finite(r$HR) & r$HR <= 0),
    missing_ids = as.character(expected_ids)[!present],
    invalid_or_failed_HR_ids = as.character(expected_ids)[!hr_ok],
    full_roster_HR_status = if (all_hr) "complete" else "unavailable_due_to_missing_failed_or_invalid_HR",
    geometric_mean_HR_primary = geo, mean_log_HR = mean_log,
    arithmetic_mean_HR_supplementary = metric(h), empirical_HR_SD = if (all_hr) sd(h) else NA_real_,
    geometric_mean_MCSE_delta_method = geo * mcse(log(h)), arithmetic_mean_MCSE = mcse(h),
    mean_squared_deviation = metric(mse_terms), squared_deviation_MCSE = mcse(mse_terms),
    squared_deviation_label = labels$squared_deviation, reference = labels$reference,
    reference_role = labels$reference_role, rejection_label = labels$rejection,
    rejection_count = reject_count, rejection_denominator_available = sum(test_ok),
    rejection_denominator_expected = n_expected,
    rejection_probability_full_roster = if (all_test) rate_available else NA_real_,
    rejection_probability_available_explicitly_conditional = rate_available,
    rejection_MCSE_available = if (any(test_ok)) sqrt(rate_available * (1 - rate_available) / sum(test_ok)) else NA_real_,
    failure_records = r[present & !fit_ok, c("replication_id", "status", "failure_reason"), drop = FALSE],
    uncertainty_note = "MCSE uses across-replication sample variation; geometric MCSE is delta-method; incomplete-denominator rates are conditional and may be selected")
}

als_baseline_summary <- function(prepared, assignments, smd_convention = NULL, config = als_config()) {
  als_validate_assignments(assignments, prepared$n, config)
  groups <- c(list(trial_treated = assignments$trial[assignments$trt == 1L],
                  trial_control = assignments$trial[assignments$trt == 0L]), assignments$external)
  summaries <- do.call(rbind, lapply(names(groups), function(g) {
    X <- prepared$raw_X[groups[[g]], , drop = FALSE]
    data.frame(group = g, covariate = config$covariates, n = nrow(X),
      mean = colMeans(X), sample_sd = apply(X, 2L, sd),
      min = apply(X, 2L, min), median = apply(X, 2L, median), max = apply(X, 2L, max),
      units = "raw_covariate_scale")
  }))
  if (is.null(smd_convention)) return(list(summaries = summaries, SMDs = NULL,
    SMD_status = "STOP_SMD_branch_pending_author_selection_of_denominator_convention"))
  als_assert(identical(smd_convention, "pooled_sample_sd"), "Unimplemented/unapproved SMD convention")
  # SMD = (mean_A - mean_B)/sqrt((sample_variance_A+sample_variance_B)/2).
  # Raw binary SEX_NUM uses this same declared sample-variance convention.
  pairs <- utils::combn(names(groups), 2L, simplify = FALSE)
  smds <- do.call(rbind, lapply(pairs, function(pair) {
    a <- summaries[summaries$group == pair[1L], ]; b <- summaries[summaries$group == pair[2L], ]
    denominator <- sqrt((a$sample_sd^2 + b$sample_sd^2) / 2)
    valid <- is.finite(denominator) & denominator > 0
    data.frame(group_A = pair[1L], group_B = pair[2L], covariate = config$covariates,
      SMD = ifelse(valid, (a$mean - b$mean) / denominator, NA_real_),
      valid = valid, denominator = denominator, convention = "mean_A_minus_mean_B_over_pooled_sample_SD")
  }))
  list(summaries = summaries, SMDs = smds,
       SMD_status = "explicit_caller_selected_convention_not_implicitly_adopted_for_ALS")
}

als_followup_summary <- function(scenario, config = als_config()) {
  groups <- c(list(trial_treated = scenario$trial[scenario$trial$trt == 1L, , drop = FALSE],
                  trial_control = scenario$trial[scenario$trial$trt == 0L, , drop = FALSE]), scenario$external)
  do.call(rbind, lapply(names(groups), function(g) {
    d <- groups[[g]]
    do.call(rbind, lapply(unique(d$outcome_origin), function(origin) {
      z <- d[d$outcome_origin == origin, , drop = FALSE]
      generated <- origin == "generated_treated"
      data.frame(scenario = scenario$scenario, group = g, outcome_origin = origin,
        censoring_origin = unique(z$censoring_origin), n = nrow(z),
        events = sum(z$status), censored = sum(z$status == 0L),
        event_proportion = mean(z$status), censoring_proportion = mean(z$status == 0L),
        admin_capped_outcomes = sum(z$administratively_truncated),
        generated_random_censorings = if (generated) sum(z$status == 0L & z$C_random < config$horizon) else NA_integer_,
        generated_admin_censorings = if (generated) sum(z$status == 0L & z$C_random >= config$horizon) else NA_integer_,
        boundary_events = sum(z$status == 1L & z$time == config$horizon),
        followup_mean_months = mean(z$time), followup_SD_months = sd(z$time),
        followup_min_months = min(z$time),
        followup_q25_months = as.numeric(quantile(z$time, 0.25, type = 7)),
        followup_median_months = median(z$time),
        followup_q75_months = as.numeric(quantile(z$time, 0.75, type = 7)),
        followup_max_months = max(z$time))
    }))
  }))
}

als_diagnostic_compatibility <- function(diagnostic, config = als_config()) {
  als_assert(identical(diagnostic$scenario, "DIAGNOSTIC"), "Diagnostic scenario required")
  followup <- als_followup_summary(diagnostic, config)
  t <- followup[followup$group == "trial_treated", ]; c <- followup[followup$group == "trial_control", ]
  list(label = "descriptive_generated_treated_vs_observed_control_compatibility_NOT_formal_type_I_error",
    event_proportion_difference = t$event_proportion - c$event_proportion,
    mean_observed_followup_difference_months = t$followup_mean_months - c$followup_mean_months,
    median_observed_followup_difference_months = t$followup_median_months - c$followup_median_months,
    treated_n = t$n, control_n = c$n, cohort_summaries = followup)
}
