# Selectively extracted definitions; no historical runner is sourced.
# See docs/PROVENANCE.csv for source identity and extraction scope.
simulate_survival_data_biased <- function(n, 
                                          source = c("trial", "external"), 
                                          ext_model_type = c("same", "different"),
                                          lambda = 0.01, 
                                          nu = 1.2, 
                                          r = 2, 
                                          delta = -0.5, 
                                          betas = NULL,
                                          total_time = 60,
                                          random_cens_rate = 0.001,
                                          cor_list = list(c(1, 2, 0.3)),
                                          # New Argument: e.g., list(var = 1, range = c(0, Inf)) 
                                          # for "all positive" bias on X1
                                          bias_config = NULL) { 
  
  source <- match.arg(source)
  ext_model_type <- match.arg(ext_model_type)
  n_vars <- if (source == "external" && ext_model_type == "different") 4 else 3
  
  # 1. Setup Correlation Matrix
  sigma_mat <- diag(n_vars)
  if (!is.null(cor_list) && !any(is.na(cor_list))) {
    for (item in cor_list) {
      idx1 <- item[1]; idx2 <- item[2]; rho  <- item[3]
      if (idx1 <= n_vars && idx2 <= n_vars) {
        sigma_mat[idx1, idx2] <- rho
        sigma_mat[idx2, idx1] <- rho
      }
    }
  }
  
  # 2. Biased Sampling via Rejection
  X_final <- matrix(NA, nrow = 0, ncol = n_vars)
  
  while(nrow(X_final) < n) {
    # Generate a batch of candidates
    batch_size <- max(n * 2, 500) 
    X_batch <- MASS::mvrnorm(n = batch_size, mu = rep(0, n_vars), Sigma = sigma_mat)
    
    if (!is.null(bias_config)) {
      # Apply the filter based on the specified range
      # bias_config$var: index of covariate (1, 2, etc.)
      # bias_config$range: c(min, max)
      keep_idx <- which(X_batch[, bias_config$var] >= bias_config$range[1] & 
                          X_batch[, bias_config$var] <= bias_config$range[2])
      X_batch <- X_batch[keep_idx, , drop = FALSE]
    }
    
    X_final <- rbind(X_final, X_batch)
  }
  
  # Trim to exactly n patients
  X_mat <- X_final[1:n, , drop = FALSE]
  colnames(X_mat) <- paste0("X", 1:n_vars)
  
  # 3. Outcome Generation (Standard Process)
  if (is.null(betas)) betas <- rep(0.5, n_vars)
  
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
  
  C_random <- if (random_cens_rate > 0) rexp(n, rate = random_cens_rate) else rep(Inf, n)
  C_admin <- rep(total_time, n)
  time_obs <- pmin(T_true, C_random, C_admin)
  status <- as.integer(T_true <= C_random & T_true <= C_admin)
  
  df_full <- cbind(
    data.frame(id = 1:n, source = source, trt = Z, time = time_obs, status = status),
    as.data.frame(X_mat)
  )
  
  return(list(data = df_full, X_mat = X_mat, trt = Z, time = time_obs, status = status))
}

