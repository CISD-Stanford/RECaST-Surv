# Selectively extracted definitions; no historical runner is sourced.
# See docs/PROVENANCE.csv for source identity and extraction scope.
simulate_survival_data <- function(n, 
                                   source = c("trial", "external"), 
                                   ext_model_type = c("same", "different"),
                                   lambda = 0.01, 
                                   nu = 1.2, 
                                   r = 2, 
                                   delta = -0.5, 
                                   betas = NULL,
                                   total_time = 60,
                                   random_cens_rate = 0.001) {
  
  source <- match.arg(source)
  ext_model_type <- match.arg(ext_model_type)
  
  if (source == "trial") {
    if (is.null(betas)) betas <- c(0.5, 0.5, 0.5)
    X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
    X_mat <- cbind(X1, X2, X3)
    colnames(X_mat) <- c("X1", "X2", "X3")
  } else { 
    if (ext_model_type == "same") {
      if (is.null(betas)) betas <- c(0.5, 0.5, 0.5)
      X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
      X_mat <- cbind(X1, X2, X3)
      colnames(X_mat) <- c("X1", "X2", "X3")
    } else {
      if (is.null(betas)) betas <- c(0.5, 0.5, 0.5, 0.5)
      X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n); X4 <- rnorm(n)
      X_mat <- cbind(X1, X2, X3, X4)
      colnames(X_mat) <- c("X1", "X2", "X3", "X4")
    }
  }
  
  if (source == "trial") {
    prob_trt <- r / (r + 1)
    Z <- rbinom(n, 1, prob_trt)
  } else {
    Z <- rep(0, n) 
  }
  
  lp <- (X_mat %*% betas) + (delta * Z)
  lp <- as.vector(lp)
  
  U <- runif(n)
  T_true <- (-log(U) / (lambda * exp(lp)))^(1/nu)
  
  if (random_cens_rate > 0) {
    C_random <- rexp(n, rate = random_cens_rate)
  } else {
    C_random <- rep(Inf, n)
  }
  C_admin <- rep(total_time, n)
  time_obs <- pmin(T_true, C_random, C_admin)
  status <- as.integer(T_true <= C_random & T_true <= C_admin)
  
  df_base <- data.frame(id = 1:n, source = source, trt = Z, time = time_obs, status = status)
  df_full <- cbind(df_base, as.data.frame(X_mat))
  
  output_list <- list(
    data = df_full, X_mat = X_mat, trt = Z, time = time_obs, status = status, t_time = T_true
  )
  return(output_list)
}

fit_external_models <- function(exter_data) {
  n_source <- length(exter_data)
  source_models <- list()
  for(ss in 1:n_source){
    source_k_data <- exter_data[[ss]]
    model <- coxph(Surv(time, status) ~ X1 + X2 + X3, data = source_k_data$data)
    source_models[[ss]] <- model
  }
  return(source_models)
}

compare_predicted_to_observed <- function(median_S_curve, time_grid, obs_time, obs_status) {
  if (is.matrix(median_S_curve)) {
    ref_S <- colMeans(median_S_curve, na.rm = TRUE)
  } else {
    ref_S <- median_S_curve
  }
  ref_H <- -log(ref_S)
  
  get_expected_hazard <- function(t) {
    approx(time_grid, ref_H, xout = t, rule = 2)$y
  }
  E_individual <- get_expected_hazard(obs_time)
  E_total <- sum(E_individual)
  O_total <- sum(obs_status)
  
  HR <- O_total / E_total
  Z_score <- (O_total - E_total) / sqrt(E_total)
  p_value <- pnorm(Z_score)
  
  return(list(O = O_total, E = E_total, HR = HR, Z = Z_score, p_value = p_value))
}

