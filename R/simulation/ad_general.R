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
                                   random_cens_rate = 0.001,
                                   # General correlation list: e.g., list(c(1, 2, 0.5), c(1, 3, 0.2))
                                   cor_list = list(c(1, 2, 0.3))) { 
  
  source <- match.arg(source)
  ext_model_type <- match.arg(ext_model_type)
  
  # Determine number of covariates
  n_vars <- if (source == "external" && ext_model_type == "different") 4 else 3
  
  # --- Define General Correlation Matrix ---
  # Initialize as Identity Matrix (Unit Variance, Zero Correlation)
  sigma_mat <- diag(n_vars)
  
  # Fill the matrix based on cor_list
  if (!is.null(cor_list) && !any(is.na(cor_list))) {
    for (item in cor_list) {
      idx1 <- item[1]
      idx2 <- item[2]
      rho  <- item[3]
      
      # Basic validation to ensure indices fit the current model dimensions
      if (idx1 <= n_vars && idx2 <= n_vars) {
        sigma_mat[idx1, idx2] <- rho
        sigma_mat[idx2, idx1] <- rho
      }
    }
  }
  
  # Safety check: Ensure the matrix is positive-definite (required for mvrnorm)
  # If a user provides an impossible correlation structure, this will catch it
  if (any(eigen(sigma_mat)$values <= 0)) {
    stop("The provided correlation matrix is not positive-definite. Please check your cor_list values.")
  }
  
  # Generate Multivariate Normal Data
  X_mat <- MASS::mvrnorm(n = n, mu = rep(0, n_vars), Sigma = sigma_mat)
  colnames(X_mat) <- paste0("X", 1:n_vars)
  
  # Set default betas
  if (is.null(betas)) betas <- rep(0.5, n_vars)
  # ----------------------------------------
  
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

