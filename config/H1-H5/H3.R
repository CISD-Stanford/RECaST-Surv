# Explicit nonproduction configuration data; author approval 2026-10-25.
# No trial identifiers, RNG assignments, runtime imports, or production task definitions.
list(
  list(
    "schema" = "H-DGM-v001",
    "status" = "NONPRODUCTION",
    "production_authorized" = FALSE,
    "scenario" = "H3",
    "variant" = "all_zero_null",
    "foundation" = "M1 Same",
    "covariates" = list(
      "names" = c("X1", "X2", "X3"),
      "distribution" = "independent standard Normal"
    ),
    "treatment_probability" = 2/3,
    "populations" = list(
      "target" = list(
        "n" = 300,
        "treatment" = "Bernoulli(2/3)",
        "beta" = c(0.2, 0.3, 0.4),
        "baseline" = list(
          "kind" = "Weibull",
          "lambda" = 0.01,
          "nu" = 1.2,
          "multiplier" = 1
        ),
        "dropout" = list(
          "kind" = "Exponential potential dropout",
          "retention_at_60" = 0.85,
          "rate" = -log(.85)/60,
          "independent_of" = c("treatment", "covariates", "event times")
        ),
        "admin_month" = 60
      ),
      "S1" = list(
        "n" = 100,
        "treatment" = "untreated",
        "beta" = c(0.2, 0.3, 0.4),
        "baseline" = list(
          "kind" = "Weibull",
          "lambda" = 0.01,
          "nu" = 1.2,
          "multiplier" = 1
        ),
        "dropout" = list(
          "kind" = "Exponential potential dropout",
          "retention_at_60" = 0.85,
          "rate" = -log(.85)/60,
          "independent_of" = c("treatment", "covariates", "event times")
        ),
        "admin_month" = 60
      ),
      "S2" = list(
        "n" = 250,
        "treatment" = "untreated",
        "beta" = c(0.15, 0.35, 0.45),
        "baseline" = list(
          "kind" = "Weibull",
          "lambda" = 0.01,
          "nu" = 1.2,
          "multiplier" = 1
        ),
        "dropout" = list(
          "kind" = "Exponential potential dropout",
          "retention_at_60" = 0.85,
          "rate" = -log(.85)/60,
          "independent_of" = c("treatment", "covariates", "event times")
        ),
        "admin_month" = 60
      ),
      "S3" = list(
        "n" = 400,
        "treatment" = "untreated",
        "beta" = c(0.22, 0.28, 0.38),
        "baseline" = list(
          "kind" = "Weibull",
          "lambda" = 0.01,
          "nu" = 1.2,
          "multiplier" = 1
        ),
        "dropout" = list(
          "kind" = "Exponential potential dropout",
          "retention_at_60" = 0.85,
          "rate" = -log(.85)/60,
          "independent_of" = c("treatment", "covariates", "event times")
        ),
        "admin_month" = 60
      )
    ),
    "effect" = list(
      "kind" = "delayed",
      "Delta" = 0,
      "switch_month" = 20,
      "pre" = 0,
      "post" = 0
    ),
    "metadata" = list(
      "changed_component" = "target treatment effect over time",
      "event_time_domain" = "all nonnegative months; no truncation at 60",
      "effect_meaning" = "generating treatment Delta; unrelated to MCMC calibration delta",
      "analysis" = "immutable original Weibull/Cauchy/GL64 kernels; not executed",
      "analysis_qualification" = "intentional working-model misspecification",
      "references" = "blocked pending later scientific authorization; no numerical values",
      "comparators" = "Cox, R-MAP, C-RMAP blocked; not implemented or activated",
      "trial_disposition" = "conditional existing-result reuse candidate; acceptance pending; no bank accessed",
      "authorized_future_new_trials" = NULL,
      "source_fits" = "future fitting to newly generated datasets; no scaling of saved fits",
      "dropout_interpretation" = "potential dropout by 60; not observed dropout or total censoring",
      "integration" = "not connected to SIM06 scenario policy, controller, workers, task or RNG allocation"
    )
  ),
  list(
    "schema" = "H-DGM-v001",
    "status" = "NONPRODUCTION",
    "production_authorized" = FALSE,
    "scenario" = "H3",
    "variant" = "delayed_alternative",
    "foundation" = "M1 Same",
    "covariates" = list(
      "names" = c("X1", "X2", "X3"),
      "distribution" = "independent standard Normal"
    ),
    "treatment_probability" = 2/3,
    "populations" = list(
      "target" = list(
        "n" = 300,
        "treatment" = "Bernoulli(2/3)",
        "beta" = c(0.2, 0.3, 0.4),
        "baseline" = list(
          "kind" = "Weibull",
          "lambda" = 0.01,
          "nu" = 1.2,
          "multiplier" = 1
        ),
        "dropout" = list(
          "kind" = "Exponential potential dropout",
          "retention_at_60" = 0.85,
          "rate" = -log(.85)/60,
          "independent_of" = c("treatment", "covariates", "event times")
        ),
        "admin_month" = 60
      ),
      "S1" = list(
        "n" = 100,
        "treatment" = "untreated",
        "beta" = c(0.2, 0.3, 0.4),
        "baseline" = list(
          "kind" = "Weibull",
          "lambda" = 0.01,
          "nu" = 1.2,
          "multiplier" = 1
        ),
        "dropout" = list(
          "kind" = "Exponential potential dropout",
          "retention_at_60" = 0.85,
          "rate" = -log(.85)/60,
          "independent_of" = c("treatment", "covariates", "event times")
        ),
        "admin_month" = 60
      ),
      "S2" = list(
        "n" = 250,
        "treatment" = "untreated",
        "beta" = c(0.15, 0.35, 0.45),
        "baseline" = list(
          "kind" = "Weibull",
          "lambda" = 0.01,
          "nu" = 1.2,
          "multiplier" = 1
        ),
        "dropout" = list(
          "kind" = "Exponential potential dropout",
          "retention_at_60" = 0.85,
          "rate" = -log(.85)/60,
          "independent_of" = c("treatment", "covariates", "event times")
        ),
        "admin_month" = 60
      ),
      "S3" = list(
        "n" = 400,
        "treatment" = "untreated",
        "beta" = c(0.22, 0.28, 0.38),
        "baseline" = list(
          "kind" = "Weibull",
          "lambda" = 0.01,
          "nu" = 1.2,
          "multiplier" = 1
        ),
        "dropout" = list(
          "kind" = "Exponential potential dropout",
          "retention_at_60" = 0.85,
          "rate" = -log(.85)/60,
          "independent_of" = c("treatment", "covariates", "event times")
        ),
        "admin_month" = 60
      )
    ),
    "effect" = list(
      "kind" = "delayed",
      "Delta" = -0.5,
      "switch_month" = 20,
      "pre" = 0,
      "post" = -0.5
    ),
    "metadata" = list(
      "changed_component" = "target treatment effect over time",
      "event_time_domain" = "all nonnegative months; no truncation at 60",
      "effect_meaning" = "generating treatment Delta; unrelated to MCMC calibration delta",
      "analysis" = "immutable original Weibull/Cauchy/GL64 kernels; not executed",
      "analysis_qualification" = "intentional working-model misspecification",
      "references" = "blocked pending later scientific authorization; no numerical values",
      "comparators" = "Cox, R-MAP, C-RMAP blocked; not implemented or activated",
      "trial_disposition" = "new trials required; not generated or allocated",
      "authorized_future_new_trials" = 1000,
      "source_fits" = "future fitting to newly generated datasets; no scaling of saved fits",
      "dropout_interpretation" = "potential dropout by 60; not observed dropout or total censoring",
      "integration" = "not connected to SIM06 scenario policy, controller, workers, task or RNG allocation"
    )
  )
)
