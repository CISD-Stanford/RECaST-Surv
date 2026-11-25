# Selectively extracted definitions; no historical runner is sourced.
# See docs/PROVENANCE.csv for source identity and extraction scope.
invert_H0 <- function(target_H, H0_table) {
  idx <- findInterval(target_H, H0_table$hazard) + 1
  idx <- pmin(pmax(idx, 1), nrow(H0_table))
  return(H0_table$time[idx])
}

build_real_survival_cache <- function(
    df,
    covariate_names,
    time_var = "STIME",
    status_var = "STATUS",
    sex_var = "SEX",
    onset_var = "ONSET"
) {
  df_std <- df
  
  if ("SEX_NUM" %in% covariate_names && !("SEX_NUM" %in% names(df_std))) {
    df_std$SEX_NUM <- ifelse(df_std[[sex_var]] == "M", 1, 0)
  }
  
  if ("ONSET_NUM" %in% covariate_names && !("ONSET_NUM" %in% names(df_std))) {
    df_std$ONSET_NUM <- as.numeric(as.factor(df_std[[onset_var]])) - 1
  }
  
  missing_vars <- setdiff(covariate_names, names(df_std))
  if (length(missing_vars) > 0) {
    stop("Missing covariates in df: ", paste(missing_vars, collapse = ", "))
  }
  
  # Standardize using FULL real dataset
  cov_means <- sapply(covariate_names, function(v) mean(df_std[[v]], na.rm = TRUE))
  cov_sds   <- sapply(covariate_names, function(v) sd(df_std[[v]], na.rm = TRUE))
  cov_sds[is.na(cov_sds) | cov_sds == 0] <- 1
  
  for (v in covariate_names) {
    df_std[[v]] <- (df_std[[v]] - cov_means[[v]]) / cov_sds[[v]]
  }
  
  # Global Cox model on standardized covariates
  cox_formula <- as.formula(
    paste0("Surv(", time_var, ", ", status_var, ") ~ ",
           paste(covariate_names, collapse = " + "))
  )
  
  global_cox <- coxph(cox_formula, data = df_std)
  H0_table   <- basehaz(global_cox, centered = FALSE)
  beta_hat   <- coef(global_cox)
  
  # Empirical exponential censoring-rate estimate
  cens_rate_empirical <- sum(df_std[[status_var]] == 0, na.rm = TRUE) /
    sum(df_std[[time_var]], na.rm = TRUE)
  
  cache <- list(
    df_std = df_std,
    covariate_names = covariate_names,
    time_var = time_var,
    status_var = status_var,
    cov_means = cov_means,
    cov_sds = cov_sds,
    global_cox = global_cox,
    H0_table = H0_table,
    beta_hat = beta_hat,
    cens_rate_empirical = cens_rate_empirical
  )
  
  return(cache)
}

sample_real_ids <- function(
    df_std,
    available_ids,
    n,
    source,
    trial_sampling = c("uniform", "younger_shorter_disdur"),
    trial_age_var = "AGE",
    trial_disdur_var = "DISDUR",
    trial_age_bias = 0.75,
    trial_disdur_bias = 0.75,
    trial_weight_cap = 20
) {
  trial_sampling <- match.arg(trial_sampling)
  
  if (source != "trial" || trial_sampling == "uniform") {
    return(sample(available_ids, size = n, replace = FALSE))
  }
  
  if (!(trial_age_var %in% names(df_std))) {
    stop("trial_age_var not found in df_std: ", trial_age_var)
  }
  if (!(trial_disdur_var %in% names(df_std))) {
    stop("trial_disdur_var not found in df_std: ", trial_disdur_var)
  }
  
  age_z    <- df_std[available_ids, trial_age_var]
  disdur_z <- df_std[available_ids, trial_disdur_var]
  
  # Because covariates are standardized in build_real_survival_cache():
  #   smaller AGE_z    = younger
  #   smaller DISDUR_z = shorter disease duration
  #
  # Positive trial_age_bias and trial_disdur_bias therefore upweight
  # younger patients and patients with shorter disease duration.
  log_w <- -trial_age_bias * age_z + trial_disdur_bias * disdur_z
  
  # Stabilize weights to avoid a few patients dominating the sample
  log_w <- log_w - max(log_w, na.rm = TRUE)
  w <- exp(log_w)
  
  med_w <- median(w, na.rm = TRUE)
  if (is.finite(med_w) && med_w > 0) {
    w <- pmin(w, trial_weight_cap * med_w)
  }
  
  w[!is.finite(w)] <- 0
  
  if (sum(w) <= 0) {
    stop("All trial sampling weights are zero or non-finite.")
  }
  
  sampled_ids <- sample(
    available_ids,
    size = n,
    replace = FALSE,
    prob = w
  )
  
  return(sampled_ids)
}

simulate_survival_data <- function(
    n,
    source = c("trial", "external"),
    ext_model_type = c("same", "different"),   # retained for compatibility
    lambda = NULL,                              # unused, retained for compatibility
    nu = NULL,                                  # unused, retained for compatibility
    r = 2,
    delta = -0.5,                               # log-HR target for treated arm
    betas = NULL,                               # unused, retained for compatibility
    total_time = 60,
    random_cens_rate = NULL,                    # if NULL, use empirical rate from cache
    cor_list = NULL,                            # unused, retained for compatibility
    cache,
    used_ids = integer(0),
    exact_randomization = TRUE,
    spike_treatment = TRUE,
    sampled_ids = NULL,
    fixed_trt = NULL,
    trial_sampling = c("uniform", "younger_shorter_disdur"),
    trial_age_var = "AGE",
    trial_disdur_var = "DISDUR",
    trial_age_bias = 0.75,
    trial_disdur_bias = 0.75,
    trial_weight_cap = 20
) {
  source <- match.arg(source)
  ext_model_type <- match.arg(ext_model_type)
  
  trial_sampling <- match.arg(trial_sampling)
  
  if (missing(cache) || is.null(cache)) {
    stop("simulate_survival_data() now requires a shared 'cache' built from real data.")
  }
  
  df_std <- cache$df_std
  covariate_names <- cache$covariate_names
  time_var <- cache$time_var
  status_var <- cache$status_var
  
  # ------------------------------------------------
  # Sample or reuse sampled_ids
  # ------------------------------------------------
  if (is.null(sampled_ids)) {
    all_ids <- seq_len(nrow(df_std))
    available_ids <- setdiff(all_ids, used_ids)
    
    if (length(available_ids) < n) {
      stop("Not enough unused real-data rows remaining to sample ", n, " observations.")
    }
    
    sampled_ids <- sample_real_ids(
      df_std = df_std,
      available_ids = available_ids,
      n = n,
      source = source,
      trial_sampling = trial_sampling,
      trial_age_var = trial_age_var,
      trial_disdur_var = trial_disdur_var,
      trial_age_bias = trial_age_bias,
      trial_disdur_bias = trial_disdur_bias,
      trial_weight_cap = trial_weight_cap
    )
  } else {
    if (length(sampled_ids) != n) {
      stop("Length of sampled_ids must equal n.")
    }
  }
  
  sub_df <- df_std[sampled_ids, , drop = FALSE]
  
  # Standardized X matrix renamed to X1, X2, ...
  X_mat <- as.matrix(sub_df[, covariate_names, drop = FALSE])
  mode(X_mat) <- "numeric"
  colnames(X_mat) <- paste0("X", seq_len(ncol(X_mat)))
  
  orig_time <- sub_df[[time_var]]
  orig_status <- sub_df[[status_var]]
  
  # ------------------------------------------------
  # Trial vs External generation
  # ------------------------------------------------
  if (source == "trial") {
    # Treatment assignment
    if (!is.null(fixed_trt)) {
      if (length(fixed_trt) != n) stop("fixed_trt must have length n.")
      Z <- as.integer(fixed_trt)
    } else {
      if (exact_randomization) {
        n_trt <- round(n * r / (r + 1))
        Z <- rep(0, n)
        Z[sample.int(n, size = n_trt, replace = FALSE)] <- 1
      } else {
        Z <- rbinom(n, 1, prob = r / (r + 1))
      }
    }
    
    # NULL SCENARIO: keep original outcomes exactly
    time_obs <- orig_time
    status_obs <- orig_status
    T_true <- rep(NA_real_, n)
    
    # ALT SCENARIO: spike in treatment arm only
    if (spike_treatment) {
      trt_idx <- which(Z == 1)
      
      if (length(trt_idx) > 0) {
        beta_use <- as.numeric(cache$beta_hat[covariate_names])
        lp_raw <- as.vector(X_mat[trt_idx, , drop = FALSE] %*% beta_use)
        
        # Cox PH spike-in using baseline hazard inversion
        U <- runif(length(trt_idx))
        target_H <- -log(U) / exp(lp_raw + delta)
        T_new <- vapply(target_H, invert_H0, numeric(1), H0_table = cache$H0_table)
        
        # Independent censoring
        cens_rate_use <- if (is.null(random_cens_rate)) {
          cache$cens_rate_empirical
        } else {
          random_cens_rate
        }
        
        if (cens_rate_use > 0) {
          C_random <- rexp(length(trt_idx), rate = cens_rate_use)
        } else {
          C_random <- rep(Inf, length(trt_idx))
        }
        
        C_admin <- rep(total_time, length(trt_idx))
        C_new <- pmin(C_random, C_admin)
        
        time_obs[trt_idx] <- pmin(T_new, C_new)
        status_obs[trt_idx] <- as.integer(T_new <= C_new)
        T_true[trt_idx] <- T_new
      }
    }
    
  } else {
    # External data always remain original observed outcomes
    Z <- rep(0, n)
    time_obs <- orig_time
    status_obs <- orig_status
    T_true <- rep(NA_real_, n)
  }
  
  df_base <- data.frame(
    id = seq_len(n),
    source = source,
    trt = Z,
    time = time_obs,
    status = status_obs
  )
  
  df_full <- cbind(df_base, as.data.frame(X_mat))
  
  output_list <- list(
    data = df_full,
    X_mat = X_mat,
    trt = Z,
    time = time_obs,
    status = status_obs,
    t_time = T_true,
    sampled_ids = sampled_ids,
    used_ids = unique(c(used_ids, sampled_ids))
  )
  
  return(output_list)
}

fit_external_models <- function(exter_data) {
  n_source <- length(exter_data)
  source_models <- list()
  
  p <- ncol(exter_data[[1]]$X_mat)
  rhs <- paste0("X", seq_len(p), collapse = " + ")
  form <- as.formula(paste("Surv(time, status) ~", rhs))
  
  for (ss in seq_len(n_source)) {
    source_k_data <- exter_data[[ss]]
    model <- coxph(form, data = source_k_data$data)
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

one_simulation <- function(trial_data, source_models, times_grid, r,
                           n_thin_draws = 1000,
                           calib_steps = 1000,
                           calib_burnin = 200,
                           calib_pred_draws = 50,
                           final_steps = 5000,
                           final_burnin = 1000) {
  
  treated_id <- trial_data$data$id[which(trial_data$data$trt == 1)]
  treated_X  <- trial_data$X_mat[treated_id, , drop = FALSE]
  
  control_id_all <- trial_data$data$id[which(trial_data$data$trt == 0)]
  n_ctrl_total   <- length(control_id_all)
  n_source       <- length(source_models)
  
  source_coefs <- lapply(source_models, coef)
  
  f_vec_all_list <- lapply(source_models, function(m) {
    ctrl_X_sub <- as.matrix(trial_data$X_mat[control_id_all, , drop = FALSE])
    as.numeric(ctrl_X_sub %*% coef(m))
  })
  
  n_calib <- 100
  
  export_objs <- c("mcmc_survival_ph", "log_post_survival_ph",
                   ".surv_int_ph", ".Q_LEG",
                   "posterior_predict_srv_prob_ph_optimized",
                   "compare_predicted_to_observed")
  
  calib_pvals <- foreach(i = 1:n_calib, .combine = c,
                         .packages = c("survival", "mvtnorm", "statmod", "matrixStats"),
                         .export = export_objs) %dopar% {
                           
                           n_train <- floor(n_ctrl_total / (r + 1))
                           train_indices <- sample(seq_along(control_id_all), n_train)
                           test_indices  <- setdiff(seq_along(control_id_all), train_indices)
                           
                           train_ids <- control_id_all[train_indices]
                           test_ids  <- control_id_all[test_indices]
                           
                           expert_medians_surv <- matrix(NA, nrow = length(times_grid), ncol = n_source)
                           expert_gammas       <- numeric(n_source)
                           
                           for (ss in 1:n_source) {
                             model_coefs <- source_coefs[[ss]]
                             f_vec_train <- f_vec_all_list[[ss]][train_indices]
                             
                             expert_recast <- mcmc_survival_ph(
                               y = trial_data$time[train_ids],
                               cens = trial_data$status[train_ids],
                               f = f_vec_train,
                               init_par = c(0, -10, 0, -2.0),
                               steps = calib_steps,
                               burnin = calib_burnin
                             )
                             
                             post_gamma_chain <- exp(expert_recast$chain[, 2])
                             expert_gammas[ss] <- median(post_gamma_chain)
                             
                             test_X <- as.matrix(trial_data$X_mat[test_ids, , drop = FALSE])
                             f_T_test <- as.numeric(test_X %*% model_coefs)
                             f_eps <- 1e-8
                             zero_idx <- abs(f_T_test) < f_eps
                             f_T_test[zero_idx] <- f_eps
                             
                             total_draws <- nrow(expert_recast$chain)
                             keep_idx <- round(seq(1, total_draws, length.out = calib_pred_draws))
                             chain_sub <- expert_recast$chain[keep_idx, , drop = FALSE]
                             
                             pred_survival <- posterior_predict_srv_prob_ph_optimized(chain_sub, f_T_test, times_grid, n_rep = 1)
                             
                             dim_p <- dim(pred_survival)[1]
                             dim_t <- dim(pred_survival)[2]
                             dim_i <- dim(pred_survival)[3]
                             dim(pred_survival) <- c(dim_p * dim_t, dim_i)
                             
                             med_vec <- matrixStats::rowMedians(pred_survival)
                             median_S_pat <- matrix(med_vec, nrow = dim_p, ncol = dim_t)
                             
                             expert_medians_surv[, ss] <- colMeans(median_S_pat)
                           }
                           
                           eps <- 1e-6
                           raw_weights <- 1 / (expert_gammas + eps)
                           norm_weights <- raw_weights / sum(raw_weights)
                           pred_median_avg <- expert_medians_surv %*% norm_weights
                           
                           stats <- compare_predicted_to_observed(
                             median_S_curve = pred_median_avg[, 1],
                             time_grid = times_grid,
                             obs_time = trial_data$data$time[test_ids],
                             obs_status = trial_data$data$status[test_ids]
                           )
                           
                           return(stats$p_value)
                         }
  
  p_val_threshold <- quantile(calib_pvals, 0.05, na.rm = TRUE)
  if (p_val_threshold > 0.05) p_val_threshold <- 0.05
  
  expert_medians_surv_final <- matrix(NA, nrow = length(times_grid), ncol = n_source)
  expert_gammas_final       <- numeric(n_source)
  
  for (ss in 1:n_source) {
    model_coefs <- source_coefs[[ss]]
    f_vec_all <- f_vec_all_list[[ss]]
    
    expert_recast <- mcmc_survival_ph(
      y = trial_data$time[control_id_all],
      cens = trial_data$status[control_id_all],
      f = f_vec_all,
      init_par = c(0, -10, 0, -2.0),
      steps = final_steps,
      burnin = final_burnin
    )
    
    expert_gammas_final[ss] <- median(exp(expert_recast$chain[, 2]))
    
    f_T_treated <- as.numeric(treated_X %*% model_coefs)
    f_eps <- 1e-8
    zero_idx <- abs(f_T_treated) < f_eps
    f_T_treated[zero_idx] <- f_eps
    
    total_draws <- nrow(expert_recast$chain)
    if (total_draws > n_thin_draws) {
      keep_idx <- round(seq(1, total_draws, length.out = n_thin_draws))
      chain_sub <- expert_recast$chain[keep_idx, , drop = FALSE]
    } else {
      chain_sub <- expert_recast$chain
    }
    
    pred_survival <- posterior_predict_srv_prob_ph_optimized(chain_sub, f_T_treated, times_grid, n_rep = 1)
    
    dim_p <- dim(pred_survival)[1]
    dim_t <- dim(pred_survival)[2]
    dim_i <- dim(pred_survival)[3]
    dim(pred_survival) <- c(dim_p * dim_t, dim_i)
    
    med_vec <- matrixStats::rowMedians(pred_survival)
    median_S_pat <- matrix(med_vec, nrow = dim_p, ncol = dim_t)
    expert_medians_surv_final[, ss] <- colMeans(median_S_pat)
  }
  
  eps <- 1e-6
  raw_weights <- 1 / (expert_gammas_final + eps)
  norm_weights <- raw_weights / sum(raw_weights)
  pred_median_final <- expert_medians_surv_final %*% norm_weights
  
  stats_final <- compare_predicted_to_observed(
    median_S_curve = pred_median_final[, 1],
    time_grid = times_grid,
    obs_time = trial_data$data$time[treated_id],
    obs_status = trial_data$data$status[treated_id]
  )
  
  analysis_data_NB <- data.frame(
    time = c(trial_data$time[treated_id], trial_data$time[control_id_all]),
    cens = c(trial_data$status[treated_id], trial_data$status[control_id_all]),
    trt  = c(rep(1, length(treated_id)), rep(0, length(control_id_all)))
  )
  
  ana_result_NB <- coxph(Surv(time, cens) ~ trt, data = analysis_data_NB)
  ana_result_NB_summ <- summary(ana_result_NB)
  
  return(list(
    final_HR = stats_final$HR,
    final_pval = stats_final$p_value,
    calibrated_threshold = as.numeric(p_val_threshold),
    nb = ana_result_NB_summ
  ))
}

