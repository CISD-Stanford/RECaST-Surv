# Additive reporting specification only; v001 numerical configuration is unchanged.
als_smd_config_v002 <- function() {
  groups <- data.frame(
    group_order = seq_len(7L),
    group = c("TRIAL_ALL", "TRIAL_TREATED", "TRIAL_CONTROL", "EXTERNAL_1",
              "EXTERNAL_2", "EXTERNAL_3", "EXTERNAL_POOLED"),
    expected_n = c(300L, 200L, 100L, 500L, 1000L, 990L, 2490L),
    legacy_group = c(NA_character_, "trial_treated", "trial_control", "external_500",
                     "external_1000", "external_990", NA_character_))
  pairs <- data.frame(
    comparison_order = seq_len(9L),
    group_A = c(rep("TRIAL_ALL", 4L), rep("TRIAL_CONTROL", 4L), "TRIAL_TREATED"),
    group_B = c("EXTERNAL_1", "EXTERNAL_2", "EXTERNAL_3", "EXTERNAL_POOLED",
                "EXTERNAL_1", "EXTERNAL_2", "EXTERNAL_3", "EXTERNAL_POOLED", "TRIAL_CONTROL"))
  pairs$comparison_id <- paste(pairs$group_A, pairs$group_B, sep = "_VS_")
  list(version = "ALS_revision_implementation_v002_SMD_extension",
       execution_scope = "synthetic_validation_only",
       covariates = c("SEX_NUM", "AGE", "BMI", "FVC", "DISDUR", "TOTAL", "SLOPE"),
       groups = groups, pairs = pairs,
       convention = "signed_mean_A_minus_mean_B_over_sqrt_unweighted_average_sample_variances",
       sample_variance_divisor = "n_minus_1", covariate_scale = "original_raw",
       sex_mean = "male_proportion", quantile_type = 7L,
       allocation_key = c("design", "replication_id"), scenario_specific = FALSE)
}
