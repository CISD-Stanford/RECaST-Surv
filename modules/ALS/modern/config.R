# Specification: author-approved ALS revised implementation v001, 2026-10-31.
# Source this file and R/als_core.R; neither file performs I/O, fitting or sampling.
als_config <- function() list(
  version = "ALS_revision_implementation_v001",
  execution_scope = "synthetic_validation_only",
  registry_candidate = NULL,
  registry_sha256 = "9b6e65cb614142d21509ffd5f9567931651284307ef9401769254f0dc1ce5563",
  units = "months", status_death = 1L, status_censored = 0L,
  covariates = c("SEX_NUM", "AGE", "BMI", "FVC", "DISDUR", "TOTAL", "SLOPE"),
  standardization = "full_registry_observed_mean_and_sample_sd",
  horizon = 120, time_grid = 1:120,
  n_registry = 2790L, n_trial = 300L, n_treated = 200L, n_control = 100L,
  n_external = c(external_500 = 500L, external_1000 = 1000L, external_990 = 990L),
  age_log_weight = -0.5, disdur_log_weight = 0.5, weight_cap_median_multiple = 20,
  censor_rate = 155 / 69897,
  scenarios = c("NULL", "ALTERNATIVE", "DIAGNOSTIC"),
  delta = c("NULL" = NA_real_, ALTERNATIVE = -0.5, DIAGNOSTIC = 0),
  nominal_alternative_reference = exp(-0.5),
  cox_ties = "efron", cox_alpha = 0.05,
  B = 100L, calibration_train = 33L, calibration_held_out = 67L,
  calibration_probability = 0.05, quantile_type = 7L, threshold_cap = 0.05,
  strict_rejection = TRUE,
  prediction_summary = "pointwise_draw_median",
  source_weight_epsilon = 1e-6, score_epsilon = 1e-8,
  interpolation = "linear_cumulative_hazard_rule_2_no_zero_anchor",
  inherited_fit = list(init_par = c(0, -10, 0, -2),
    calibration_steps = 1000L, calibration_burnin = 200L, calibration_draws = 50L,
    final_steps = 5000L, final_burnin = 1000L, final_draws = 1000L, n_rep = 1L),
  # No registry loader, production seed, roster, runner or job submission exists.
  real_data_execution_authorized = FALSE
)
