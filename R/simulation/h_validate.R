# H development v001: independent approved-parameter and input checks.
# No simulation, random sampling, fitting, side effects, or inherited imports.
H_require <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  invisible(TRUE)
}

H_number <- function(x, name, lower = -Inf, upper = Inf) {
  H_require(is.numeric(x) && is.null(dim(x)) && length(x) > 0L &&
              all(is.finite(x)) && all(x >= lower & x <= upper),
            paste0(name, " must be a nonempty finite numeric vector in range"))
  invisible(TRUE)
}

H_scalar <- function(x, name, lower = -Inf, upper = Inf) {
  H_number(x, name, lower, upper)
  H_require(length(x) == 1L, paste0(name, " must be scalar"))
  invisible(TRUE)
}

H_keys <- function(x, expected, name) {
  H_require(is.list(x) && identical(sort(names(x)), sort(expected)),
            paste0(name, " has missing, extra, or duplicate fields"))
}

H_validate_config <- function(cfg) {
  H_keys(cfg, c("schema", "status", "production_authorized", "scenario", "variant",
               "foundation", "covariates", "treatment_probability", "populations",
               "effect", "metadata"), "config")
  H_require(identical(cfg$schema, "H-DGM-v001") &&
              identical(cfg$status, "NONPRODUCTION") &&
              identical(cfg$production_authorized, FALSE) &&
              identical(cfg$foundation, "M1 Same"), "Invalid scope or foundation")
  H_require(is.character(cfg$scenario) && length(cfg$scenario) == 1L &&
              cfg$scenario %in% c("M1", paste0("H", 1:5)), "Invalid scenario")
  H_require(identical(cfg$covariates,
    list(names = c("X1", "X2", "X3"), distribution = "independent standard Normal")),
    "Covariate law differs from approved M1")
  H_require(identical(cfg$treatment_probability, 2/3), "Treatment probability must be 2/3")
  H_keys(cfg$populations, c("target", "S1", "S2", "S3"), "populations")
  sc <- cfg$scenario
  H_keys(cfg$effect, c("kind", "Delta", "switch_month", "pre", "post"), "effect")
  ef <- cfg$effect
  H_scalar(ef$Delta, "generating Delta")
  H_require(ef$Delta %in% c(0, -0.5), "Unapproved generating Delta")
  if (sc == "H3") {
    H_require(identical(ef$kind, "delayed") && identical(ef$switch_month, 20) &&
      identical(ef$pre, 0) && identical(ef$post, ef$Delta), "Invalid H3 effect")
  } else {
    H_require(identical(ef$kind, "constant") && is.null(ef$switch_month) &&
      identical(ef$pre, ef$Delta) && identical(ef$post, ef$Delta), "Invalid PH effect")
  }
  expected_beta <- list(target = c(.20, .30, .40), S1 = c(.20, .30, .40),
                        S2 = c(.15, .35, .45), S3 = c(.22, .28, .38))
  # Explicit approved H4 table, independent of any configuration constructor.
  if (sc == "H4") expected_beta <- list(target = c(.10, .15, .20),
    S1 = c(.10, .15, .20), S2 = c(.075, .175, .225), S3 = c(.11, .14, .19))
  sizes <- c(target = 300, S1 = 100, S2 = 250, S3 = 400)
  for (pop in names(sizes)) {
    p <- cfg$populations[[pop]]
    H_keys(p, c("n", "treatment", "beta", "baseline", "dropout", "admin_month"), pop)
    H_require(identical(p$n, unname(sizes[pop])) &&
      identical(p$treatment, if (pop == "target") "Bernoulli(2/3)" else "untreated") &&
      identical(p$beta, expected_beta[[pop]]) && identical(p$admin_month, 60),
      paste0("Unapproved population parameters: ", pop))
    b <- p$baseline
    if (sc == "H2") {
      H_keys(b, c("kind", "switch_month", "a", "b", "continuation"), "baseline")
      H_require(identical(b$kind, "piecewise_exponential") &&
        identical(b$switch_month, 30) && identical(b$a, (.01 * 60^1.2)/90) &&
        identical(b$b, 2*((.01 * 60^1.2)/90)) &&
        identical(b$continuation, "post-switch slope for all t > 30"), "Invalid H2 baseline")
    } else {
      H_keys(b, c("kind", "lambda", "nu", "multiplier"), "baseline")
      H_require(identical(b$kind, "Weibull") && identical(b$lambda, .01) &&
        identical(b$nu, 1.2), "Invalid Weibull baseline")
      H_scalar(b$multiplier, "baseline multiplier")
      H_require(if (sc == "H1" && pop != "target") b$multiplier %in% c(.5, 1, 2)
                else identical(b$multiplier, 1), "Unapproved baseline multiplier")
    }
    d <- p$dropout
    retention <- if (sc == "H5" && pop == "target") .60 else .85
    H_keys(d, c("kind", "retention_at_60", "rate", "independent_of"), "dropout")
    H_require(identical(d$kind, "Exponential potential dropout") &&
      identical(d$retention_at_60, retention) && identical(d$rate, -log(retention)/60) &&
      identical(d$independent_of, c("treatment", "covariates", "event times")),
      paste0("Invalid dropout mechanism: ", pop))
  }
  m <- if (sc == "H1") cfg$populations$S1$baseline$multiplier else 1
  if (sc == "H1") H_require(all(vapply(cfg$populations[c("S1", "S2", "S3")],
    function(p) identical(p$baseline$multiplier, m), logical(1))), "H1 needs one common source multiplier")
  variant <- if (sc == "H1") paste0("m", m, "_Delta", ef$Delta) else
    if (sc == "H3") if (ef$Delta == 0) "all_zero_null" else "delayed_alternative" else
      paste0("Delta", ef$Delta)
  H_require(identical(cfg$variant, variant), "Variant disagrees with DGM")
  reuse <- (sc == "H1" && m == 1) || (sc == "H3" && ef$Delta == 0)
  changed <- c(M1 = "none", H1 = "source baseline multiplier", H2 = "common event baseline",
               H3 = "target treatment effect over time", H4 = "all prognostic coefficients",
               H5 = "target potential dropout")
  expected_metadata <- list(
    changed_component = unname(changed[sc]),
    event_time_domain = "all nonnegative months; no truncation at 60",
    effect_meaning = "generating treatment Delta; unrelated to MCMC calibration delta",
    analysis = "immutable original Weibull/Cauchy/GL64 kernels; not executed",
    analysis_qualification = if (sc %in% c("H2", "H3")) "intentional working-model misspecification" else "inherited working analysis",
    references = "blocked pending later scientific authorization; no numerical values",
    comparators = "Cox, R-MAP, C-RMAP blocked; not implemented or activated",
    trial_disposition = if (sc == "M1") "foundation formula only; no trial bank accessed" else
      if (reuse) "conditional existing-result reuse candidate; acceptance pending; no bank accessed" else "new trials required; not generated or allocated",
    authorized_future_new_trials = if (sc == "H3" && ef$Delta == -.5) 1000 else NULL,
    source_fits = "future fitting to newly generated datasets; no scaling of saved fits",
    dropout_interpretation = "potential dropout by 60; not observed dropout or total censoring",
    integration = "not connected to SIM06 scenario policy, controller, workers, task or RNG allocation")
  H_require(identical(cfg$metadata, expected_metadata), "Scientific metadata disagrees with DGM or scope")
  invisible(TRUE)
}

H_population <- function(cfg, population) {
  H_validate_config(cfg)
  H_require(is.character(population) && length(population) == 1L &&
    population %in% names(cfg$populations), "Population must be target, S1, S2 or S3")
  cfg$populations[[population]]
}

H_X <- function(X) {
  if (is.numeric(X) && is.null(dim(X)) && length(X) == 3L) X <- matrix(X, nrow = 1L)
  H_require(is.matrix(X) && is.numeric(X) && nrow(X) > 0L && ncol(X) == 3L &&
    all(is.finite(X)), "X must have three finite numeric columns")
  H_require(is.null(colnames(X)) || identical(colnames(X), c("X1", "X2", "X3")),
    "X column order must be X1, X2, X3")
  X
}

H_align <- function(...) {
  x <- list(...); n <- max(lengths(x))
  H_require(all(lengths(x) %in% c(1L, n)), "Incompatible vector lengths; partial recycling forbidden")
  lapply(x, rep_len, length.out = n)
}

H_arm <- function(A, population) {
  H_number(A, "A", 0, 1)
  H_require(all(A %in% c(0, 1)) && (population == "target" || all(A == 0)),
    "A must be binary; external populations are untreated")
  invisible(TRUE)
}
