# Exact selected definitions from audited Code/MCMC_PH.R.
# No historical runner, library calls or execution block copied.

.Q_LEG  <- local({
  q <- statmod::gauss.quad(n = 64, kind = "legendre") 
  a <- -pi/2; b <- pi/2
  t_nodes <- ( (b - a) * (q$nodes + 1) / 2 ) + a
  t_weights <- q$weights * (b - a) / 2
  list(t = t_nodes, wt = t_weights)
})

.surv_int_ph <- function(y, cens, f, delta, gamma, nu, alpha) {
  # y: vector of observed TIMES (not log-times)
  # cens: vector of status (1=Event, 0=Censored)
  # f: vector of structural predictions
  
  t_nodes <- .Q_LEG$t
  wt      <- .Q_LEG$wt / pi
  
  # Map quadrature nodes to beta space
  tan_t <- tan(t_nodes)
  beta_nodes <- gamma * tan_t + delta
  
  # Calculate linear predictor matrix (nodes x observations)
  # eta_ij = beta_i * f_j
  Eta <- outer(beta_nodes, f) + alpha
  
  # Prepare y matrix
  Y_mat <- matrix(rep(y, each = length(t_nodes)), nrow = length(t_nodes))
  
  # Precompute common terms
  # Hazard scale: lambda = exp(eta)
  Exp_Eta <- exp(Eta)
  
  # Cumulative Hazard Base: y^nu
  # We use log-space for stability to calculate y^nu
  # y^nu = exp(nu * log(y))
  Y_pow_nu <- exp(nu * log(Y_mat))
  
  # Cumulative Hazard: H(t) = y^nu * exp(eta)
  H_t <- Y_pow_nu * Exp_Eta
  
  # Survival Function: S(t) = exp(-H(t))
  Surv_val <- exp(-H_t)
  
  # Density (for events): f(t) = h(t) * S(t)
  # h(t) = nu * y^(nu-1) * exp(eta)
  #      = (nu/y) * (y^nu * exp(eta)) = (nu/y) * H_t
  # This avoids recomputing powers, but y must be > 0
  # Add small epsilon to y to prevent div/0 if y=0 exists (unlikely in continuous)
  Haz_val <- (nu / (Y_mat + 1e-10)) * H_t
  Dens_val <- Haz_val * Surv_val
  
  # Select Likelihood based on censoring
  Cens_mat <- matrix(rep(cens, each = length(t_nodes)), nrow = length(t_nodes))
  Likelihood <- ifelse(Cens_mat == 1, Dens_val, Surv_val)
  
  # Integration: Weighted sum
  drop(as.vector(wt %*% Likelihood))
}

log_post_survival_ph <- function(par, y, cens, f) {
  delta     <- par[1]
  log_gamma <- par[2]
  log_nu    <- par[3] # Weibull Shape
  alpha     <- par[4]
  
  gamma <- exp(log_gamma)
  nu    <- exp(log_nu)
  n_T   <- length(y)
  
  integ <- .surv_int_ph(y, cens, f, delta, gamma, nu, alpha)
  
  # Stability check
  if (any(!is.finite(integ)) || any(integ <= 1e-300)) return(-Inf)
  
  ll <- sum(log(integ))
  
  # Priors (scaled by 1/n_T)
  # nu ~ LogNormal(0, 1) -> log_nu ~ Normal(0, 1)
  lp <- ll +
    dnorm(delta,     1, 20, log = TRUE) / n_T +
    dnorm(log_gamma, 0,  3, log = TRUE) / n_T + 
    dnorm(log_nu,    0,  1, log = TRUE) / n_T +
    dnorm(alpha,     -4, 10, log = TRUE) / n_T    # <--- Prior for Alpha
  
  return(lp)
}

mcmc_survival_ph <- function(y, cens, f, init_par, steps = 5000, burnin = 1000) {
  par    <- as.numeric(init_par)
  n_par  <- length(par)
  chain  <- matrix(NA_real_, nrow = steps, ncol = n_par)
  colnames(chain) <- c("delta", "log_gamma", "log_nu", "alpha")
  
  pcov   <- diag(n_par)
  pscale <- 1e-4
  accept <- 0L
  
  # Initial check
  log_dens_prev <- log_post_survival_ph(par, y, cens, f)
  if (!is.finite(log_dens_prev)) stop("Bad initial parameters (log posterior = -Inf).")
  
  chain[1, ] <- par
  
  for (ttt in 2:steps) {
    # Proposal
    prop <- as.numeric(mvtnorm::rmvnorm(1, mean = par, sigma = pcov * pscale))
    log_dens <- log_post_survival_ph(prop, y, cens, f)
    
    # Early chain robustness
    if (ttt < burnin) {
      tries <- 0L
      while (!is.finite(log_dens) && tries < 20) {
        prop <- as.numeric(mvtnorm::rmvnorm(1, mean = par, sigma = pcov * pscale))
        log_dens <- log_post_survival_ph(prop, y, cens, f)
        tries <- tries + 1L
      }
    }
    
    # Metropolis Step
    if (is.finite(log_dens) && (log_dens - log_dens_prev) > log(runif(1))) {
      par <- prop
      log_dens_prev <- log_dens
      accept <- accept + 1L
    }
    chain[ttt, ] <- par
    
    # Adaptation
    if (ttt < burnin) {
      if (ttt == 100) pscale <- 0.1
      tmp <- if (ttt >= 200) chain[(ttt-199):ttt, , drop = FALSE] else NULL
      if (!is.null(tmp)) {
        tmpu <- unique(tmp)
        if (nrow(tmpu) >= 4) pcov <- stats::cov(tmpu) + diag(n_par)*1e-6
      }
      if (ttt %% 50 == 0) {
        rate <- accept / 50
        accept <- 0L
        if (rate < 0.2) pscale <- pscale * 0.8
        if (rate > 0.3) pscale <- pscale * 1.2
      }
    }
  }
  list(chain = chain[burnin:steps, , drop = FALSE], 
       final_accept_rate = accept / (steps - burnin))
}

posterior_predict_srv_prob_ph_optimized <- function(chain, f_treated, time_grid, n_rep = 10) {
  
  # 1. Setup Dimensions
  n_draws <- nrow(chain)
  n_pats  <- length(f_treated)
  n_times <- length(time_grid)
  total_iters <- n_draws * n_rep
  
  # 2. Expand Parameters (Vectorization)
  # We replicate the rows of the chain to match total_iters
  # idx = c(1, 1, ..., 2, 2, ...)
  idx <- rep(seq_len(n_draws), each = n_rep)
  
  delta_all <- chain[idx, 1]
  gamma_all <- exp(chain[idx, 2])
  nu_all    <- exp(chain[idx, 3])
  alpha_all <- chain[idx, 4]
  
  # 3. Generate Random Effects (All at once)
  # We generate n_pats * total_iters random numbers in one go.
  # The location/scale vectors are repeated 'each=n_pats' to match the 
  # column-major filling of the matrix (Patient 1..P for Iter 1, then Iter 2...)
  beta_vals <- rcauchy(n_pats * total_iters, 
                       location = rep(delta_all, each = n_pats), 
                       scale    = rep(gamma_all, each = n_pats))
  
  # Reshape into Matrix [n_pats x total_iters]
  beta_mat <- matrix(beta_vals, nrow = n_pats, ncol = total_iters)
  
  # 4. Calculate Hazard Scale
  # LP = alpha + beta * f
  # Note: 'beta_mat * f_treated' recycles f_treated column-wise (correct)
  # We calculate the linear predictor for every patient/iteration simultaneously
  # t(t(...) + alpha_all) adds alpha_all[j] to the j-th column
  lp_mat <- t(t(beta_mat * f_treated) + alpha_all)
  
  # Hazard Scale = exp(LP) -> Matrix [n_pats x total_iters]
  hazard_scale_mat <- exp(lp_mat)
  
  # 5. Compute Probabilities over Time Grid
  # We loop over time points (usually small, < 100) instead of iterations (thousands)
  out_probs <- array(NA_real_, dim = c(n_pats, n_times, total_iters))
  
  for (j in seq_len(n_times)) {
    t_val <- time_grid[j]
    
    # Calculate t^nu for all iterations (Vector length total_iters)
    # This handles the varying 'nu' per iteration
    time_term <- t_val ^ nu_all
    
    # Cumulative Hazard = HazardScale * time_term
    # Broadcast multiply: Column j of hazard_scale * Element j of time_term
    cum_haz_slice <- t(t(hazard_scale_mat) * time_term)
    
    # Survival Probability S(t) = exp(-H(t))
    out_probs[, j, ] <- exp(-cum_haz_slice)
  }
  
  return(out_probs)
}
