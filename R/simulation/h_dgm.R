# Pure arithmetic only. Requires validate.R in the same environment.
H_linear_predictor <- function(cfg, population, X) {
  p <- H_population(cfg, population)
  lp <- as.vector(H_X(X) %*% p$beta)
  H_number(lp, "linear predictor")
  lp
}

H_baseline_cumulative <- function(cfg, population, t) {
  b <- H_population(cfg, population)$baseline
  H_number(t, "t", 0)
  if (b$kind == "Weibull") b$multiplier*b$lambda*t^b$nu else
    b$a*pmin(t, b$switch_month) + b$b*pmax(t-b$switch_month, 0)
}

H_baseline_hazard <- function(cfg, population, t) {
  b <- H_population(cfg, population)$baseline
  H_number(t, "t", 0)
  if (b$kind == "Weibull") b$multiplier*b$lambda*b$nu*t^(b$nu-1) else
    ifelse(t <= b$switch_month, b$a, b$b)
}

H_inverse_baseline <- function(cfg, population, q) {
  b <- H_population(cfg, population)$baseline
  H_number(q, "cumulative hazard input", 0)
  if (b$kind == "Weibull") (q/(b$multiplier*b$lambda))^(1/b$nu) else {
    boundary <- b$a*b$switch_month
    ifelse(q <= boundary, q/b$a, b$switch_month+(q-boundary)/b$b)
  }
}

H_treatment_log_hr <- function(cfg, t) {
  H_validate_config(cfg); H_number(t, "t", 0)
  ef <- cfg$effect
  if (ef$kind == "constant") rep(ef$Delta, length(t)) else
    ifelse(t <= ef$switch_month, ef$pre, ef$post)
}

H_hazard <- function(cfg, population, t, X, A) {
  H_arm(A, population)
  z <- H_align(t = t, lp = H_linear_predictor(cfg, population, X), A = A)
  H_baseline_hazard(cfg, population, z$t) * exp(z$lp+z$A*H_treatment_log_hr(cfg, z$t))
}

H_cumulative_hazard <- function(cfg, population, t, X, A) {
  H_arm(A, population)
  z <- H_align(t = t, lp = H_linear_predictor(cfg, population, X), A = A)
  ef <- cfg$effect
  full <- H_baseline_cumulative(cfg, population, z$t)
  if (ef$kind == "constant") full*exp(z$lp+z$A*ef$Delta) else {
    early <- H_baseline_cumulative(cfg, population, pmin(z$t, ef$switch_month))
    exp(z$lp)*(exp(z$A*ef$pre)*early + exp(z$A*ef$post)*(full-early))
  }
}

H_inverse_hazard <- function(cfg, population, q, X, A) {
  H_number(q, "cumulative hazard input", 0); H_arm(A, population)
  z <- H_align(q = q, lp = H_linear_predictor(cfg, population, X), A = A)
  ef <- cfg$effect
  if (ef$kind == "constant") {
    value <- H_inverse_baseline(cfg, population, z$q/exp(z$lp+z$A*ef$Delta))
  } else {
    base_switch <- H_baseline_cumulative(cfg, population, ef$switch_month)
    early_scale <- exp(z$A*ef$pre); late_scale <- exp(z$A*ef$post)
    w <- z$q/exp(z$lp); boundary <- early_scale*base_switch
    q_base <- ifelse(w <= boundary, w/early_scale, base_switch+(w-boundary)/late_scale)
    value <- H_inverse_baseline(cfg, population, q_base)
  }
  H_number(value, "inverse event time", 0)
  value
}

H_event_time_from_uniform <- function(cfg, population, u, X, A) {
  H_number(u, "event Uniform input", 0, 1)
  H_require(all(u > 0 & u < 1), "event Uniform inputs must lie strictly inside (0,1)")
  H_inverse_hazard(cfg, population, -log(u), X, A)
}

H_dropout_rate <- function(cfg, population) H_population(cfg, population)$dropout$rate

H_dropout_probability <- function(cfg, population, t) {
  H_number(t, "t", 0)
  -expm1(-H_dropout_rate(cfg, population)*t)
}

# q is a supplied nonnegative unit-rate Exponential latent numerical input.
# No treatment, X, or event-time argument can influence potential dropout.
H_potential_dropout_time <- function(cfg, population, q) {
  H_number(q, "dropout Exponential input", 0)
  q/H_dropout_rate(cfg, population)
}

H_dropout_time_from_uniform <- function(cfg, population, u) {
  H_number(u, "dropout Uniform input", 0, 1)
  H_require(all(u > 0 & u < 1), "dropout Uniform inputs must lie strictly inside (0,1)")
  H_potential_dropout_time(cfg, population, -log(u))
}

# Inherited event-at-tie convention; no patient/trial identities are assigned.
H_observation <- function(cfg, population, event_time, dropout_time) {
  p <- H_population(cfg, population)
  H_number(event_time, "event time", 0); H_number(dropout_time, "dropout time", 0)
  z <- H_align(event = event_time, dropout = dropout_time)
  list(time = pmin(z$event, z$dropout, p$admin_month),
       status = as.integer(z$event <= z$dropout & z$event <= p$admin_month))
}
