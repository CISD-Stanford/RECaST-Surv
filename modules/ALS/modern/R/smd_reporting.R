# Requires unchanged v001 config/core definitions and v002 config.R.
# No data I/O, model calls, RNG, scenario outcome generation or job execution.
als_smd_key_v002 <- function(...) paste(..., sep = "\034")

als_smd_id_v002 <- function(replication_id, design) {
  als_assert(length(replication_id) == 1L && !is.na(replication_id) &&
    nzchar(as.character(replication_id)) && !grepl("\034", as.character(replication_id), fixed = TRUE),
    "A nonmissing, nonempty allocation replication_id is required")
  als_assert(length(design) == 1L && design %in% c("random", "biased"), "Unknown sampling design")
  invisible(TRUE)
}

als_smd_groups_v002 <- function(prepared, assignments, context) {
  als_require_synthetic(context)
  als_validate_assignments(assignments, prepared$n, als_config())
  als_assert(is.matrix(prepared$raw_X) && is.numeric(prepared$raw_X) &&
    nrow(prepared$raw_X) == prepared$n && !anyDuplicated(colnames(prepared$raw_X)) &&
    all(als_smd_config_v002()$covariates %in% colnames(prepared$raw_X)),
    "Aligned original raw covariates are required")
  groups <- list(TRIAL_ALL = assignments$trial,
    TRIAL_TREATED = assignments$trial[assignments$trt == 1L],
    TRIAL_CONTROL = assignments$trial[assignments$trt == 0L],
    EXTERNAL_1 = assignments$external[[1L]], EXTERNAL_2 = assignments$external[[2L]],
    EXTERNAL_3 = assignments$external[[3L]],
    EXTERNAL_POOLED = unlist(assignments$external, use.names = FALSE))
  cfg <- als_smd_config_v002()
  als_assert(identical(names(groups), cfg$groups$group) &&
    identical(as.integer(lengths(groups)), cfg$groups$expected_n), "Seven-group size/order mismatch")
  groups
}

# Incomplete/nonfinite values are counted and make that group/covariate summary
# unavailable. There is no complete-case deletion, imputation or SD replacement.
als_smd_describe_v002 <- function(x, binary = FALSE) {
  n <- length(x)
  n_missing <- sum(is.na(x))
  n_nonfinite <- sum(!is.na(x) & !is.finite(x))
  n_invalid_binary <- if (binary) sum(is.finite(x) & !x %in% c(0, 1)) else 0L
  mean_valid <- n > 0L && n_missing == 0L && n_nonfinite == 0L && n_invalid_binary == 0L
  m <- if (mean_valid) mean(x) else NA_real_
  mean_valid <- mean_valid && is.finite(m)
  variance <- if (mean_valid && n >= 2L) stats::var(x) else NA_real_
  sd_valid <- is.finite(variance) && variance >= 0
  reason <- if (n_missing > 0L) "missing_values" else if (n_nonfinite > 0L) "nonfinite_values" else
    if (n_invalid_binary > 0L) "invalid_SEX_NUM_values" else if (!mean_valid) "undefined_mean" else
    if (!sd_valid) "undefined_sample_variance" else "ok"
  list(n = n, n_observed_finite = sum(is.finite(x)), n_missing = n_missing,
    n_nonfinite = n_nonfinite, n_invalid_binary = n_invalid_binary,
    mean = if (mean_valid) m else NA_real_, sample_variance = variance,
    sample_sd = if (sd_valid) sqrt(variance) else NA_real_,
    male_proportion = if (binary && mean_valid) m else NA_real_,
    min = if (mean_valid) min(x) else NA_real_, median = if (mean_valid) median(x) else NA_real_,
    max = if (mean_valid) max(x) else NA_real_, mean_valid = mean_valid,
    sample_sd_valid = sd_valid, valid = mean_valid && sd_valid, validity_reason = reason)
}

# Call once per replication/design, without a scenario argument. The optional
# state interface below provides idempotent reuse by all three paired scenarios.
als_baseline_summary_v002 <- function(prepared, assignments, replication_id, context) {
  als_require_synthetic(context)
  als_smd_id_v002(replication_id, assignments$design)
  cfg <- als_smd_config_v002()
  groups <- als_smd_groups_v002(prepared, assignments, context)
  summaries <- do.call(rbind, lapply(seq_len(nrow(cfg$groups)), function(g) {
    do.call(rbind, lapply(seq_along(cfg$covariates), function(j) {
      covariate <- cfg$covariates[j]
      stats <- als_smd_describe_v002(prepared$raw_X[groups[[g]], covariate], covariate == "SEX_NUM")
      data.frame(replication_id = as.character(replication_id), design = assignments$design,
        group_order = g, group = cfg$groups$group[g], legacy_group = cfg$groups$legacy_group[g],
        covariate_order = j, covariate = covariate, as.data.frame(stats),
        units = "raw_covariate_scale", row.names = NULL)
    }))
  }))
  smds <- do.call(rbind, lapply(seq_len(nrow(cfg$pairs)), function(k) {
    pair <- cfg$pairs[k, ]
    a <- summaries[summaries$group == pair$group_A, , drop = FALSE]
    b <- summaries[summaries$group == pair$group_B, , drop = FALSE]
    denominator <- sqrt((a$sample_variance + b$sample_variance) / 2)
    valid <- a$valid & b$valid & is.finite(denominator) & denominator > 0
    signed <- rep(NA_real_, length(cfg$covariates))
    signed[valid] <- (a$mean[valid] - b$mean[valid]) / denominator[valid]
    valid <- valid & is.finite(signed)
    signed[!valid] <- NA_real_
    reasons <- vapply(seq_along(valid), function(j) {
      if (!a$valid[j]) return(paste0("group_A_", a$validity_reason[j]))
      if (!b$valid[j]) return(paste0("group_B_", b$validity_reason[j]))
      if (!is.finite(denominator[j])) return("undefined_pooled_sample_SD")
      if (denominator[j] == 0) return("zero_pooled_sample_SD")
      if (!valid[j]) return("nonfinite_SMD")
      "ok"
    }, character(1))
    data.frame(replication_id = as.character(replication_id), design = assignments$design,
      comparison_order = k, comparison_id = pair$comparison_id,
      group_A = pair$group_A, group_B = pair$group_B,
      group_order_A = match(pair$group_A, cfg$groups$group),
      group_order_B = match(pair$group_B, cfg$groups$group),
      covariate_order = seq_along(cfg$covariates), covariate = cfg$covariates,
      n_A = a$n, n_B = b$n, mean_A = a$mean, mean_B = b$mean,
      sample_sd_A = a$sample_sd, sample_sd_B = b$sample_sd,
      sample_variance_A = a$sample_variance, sample_variance_B = b$sample_variance,
      male_proportion_A = a$male_proportion, male_proportion_B = b$male_proportion,
      n_missing_A = a$n_missing, n_missing_B = b$n_missing,
      n_nonfinite_A = a$n_nonfinite, n_nonfinite_B = b$n_nonfinite,
      n_invalid_binary_A = a$n_invalid_binary, n_invalid_binary_B = b$n_invalid_binary,
      group_A_valid = a$valid, group_B_valid = b$valid,
      pooled_sample_SD = denominator, signed_SMD = signed, absolute_SMD = abs(signed),
      valid = valid, validity_reason = reasons, convention = cfg$convention, units = "raw_covariate_scale",
      row.names = NULL)
  }))
  rownames(summaries) <- rownames(smds) <- NULL
  list(version = cfg$version, replication_id = as.character(replication_id), design = assignments$design,
       groups = groups, summaries = summaries, SMDs = smds,
       scenario_specific = FALSE, SMD_status = "author_approved_pooled_sample_SD_v002")
}

als_smd_new_state_v002 <- function(context) {
  als_require_synthetic(context)
  e <- new.env(parent = emptyenv())
  e$entries <- list(); e$calculation_count <- 0L
  e
}

als_smd_once_v002 <- function(state, prepared, assignments, replication_id, context) {
  als_require_synthetic(context)
  als_smd_id_v002(replication_id, assignments$design)
  key <- als_smd_key_v002(assignments$design, as.character(replication_id))
  # Outcomes, generated draws and scenario labels intentionally do not enter.
  snapshot <- list(n = prepared$n, raw_X = prepared$raw_X[, als_smd_config_v002()$covariates, drop = FALSE],
                   trial = assignments$trial, trt = assignments$trt, external = assignments$external)
  if (!is.null(state$entries[[key]])) {
    als_assert(identical(snapshot, state$entries[[key]]$snapshot),
               "Same replication/design has changed membership or raw covariates")
    return(state$entries[[key]]$bundle)
  }
  bundle <- als_baseline_summary_v002(prepared, assignments, replication_id, context)
  state$entries[[key]] <- list(snapshot = snapshot, bundle = bundle)
  state$calculation_count <- state$calculation_count + 1L
  bundle
}

als_smd_records_v002 <- function(state) {
  if (length(state$entries)) return(do.call(rbind, lapply(state$entries, function(e) e$bundle$SMDs)))
  data.frame(replication_id = character(), design = character(), comparison_id = character(),
    group_A = character(), group_B = character(), covariate = character(),
    signed_SMD = numeric(), absolute_SMD = numeric(), valid = logical())
}

# Quantiles use valid allocation records, with all unavailable/missing counts
# exposed against a separately supplied design/replication roster. No scenario
# row is accepted, deduplicated or counted as an independent allocation.
als_smd_across_replications_v002 <- function(records, expected_allocations, context) {
  als_require_synthetic(context)
  cfg <- als_smd_config_v002()
  als_assert(is.data.frame(expected_allocations) &&
    all(c("replication_id", "design") %in% names(expected_allocations)) && nrow(expected_allocations) > 0L &&
    !"scenario" %in% names(expected_allocations), "An allocation-only expected roster is required")
  for (i in seq_len(nrow(expected_allocations)))
    als_smd_id_v002(expected_allocations$replication_id[i], expected_allocations$design[i])
  roster_keys <- with(expected_allocations, als_smd_key_v002(design, as.character(replication_id)))
  als_assert(!anyDuplicated(roster_keys), "Duplicate replication/design in expected roster")
  required <- c("replication_id", "design", "comparison_id", "group_A", "group_B", "covariate",
                "signed_SMD", "absolute_SMD", "valid")
  als_assert(is.data.frame(records) && all(required %in% names(records)) &&
    !"scenario" %in% names(records), "SMD records must be allocation-only; scenario rows are forbidden")
  if (nrow(records)) {
    allocation_keys <- with(records, als_smd_key_v002(design, as.character(replication_id)))
    pair_index <- match(records$comparison_id, cfg$pairs$comparison_id)
    als_assert(!anyNA(pair_index) && all(allocation_keys %in% roster_keys) &&
      all(records$covariate %in% cfg$covariates) &&
      all(records$group_A == cfg$pairs$group_A[pair_index]) &&
      all(records$group_B == cfg$pairs$group_B[pair_index]), "Unexpected allocation/comparison/covariate record")
    als_assert(is.logical(records$valid) && !anyNA(records$valid), "Explicit validity flags are required")
    valid <- records$valid
    als_assert(all(is.finite(records$signed_SMD[valid])) &&
      all(is.finite(records$absolute_SMD[valid])) &&
      all(records$absolute_SMD[valid] == abs(records$signed_SMD[valid])),
      "Valid signed/absolute SMD values must be finite and consistent")
  }
  row_keys <- with(records, als_smd_key_v002(design, as.character(replication_id), comparison_id, covariate))
  als_assert(!anyDuplicated(row_keys), "Duplicate allocation/comparison/covariate; paired scenarios must not be stacked")
  designs <- c("random", "biased")[c("random", "biased") %in% expected_allocations$design]
  result <- do.call(rbind, lapply(designs, function(design) {
    ids <- as.character(expected_allocations$replication_id[expected_allocations$design == design])
    do.call(rbind, lapply(seq_len(nrow(cfg$pairs)), function(k) {
      pair <- cfg$pairs[k, ]
      do.call(rbind, lapply(seq_along(cfg$covariates), function(j) {
        covariate <- cfg$covariates[j]
        index <- match(als_smd_key_v002(design, ids, pair$comparison_id, covariate), row_keys)
        present <- !is.na(index)
        good <- rep(FALSE, length(ids))
        good[present] <- records$valid[index[present]]
        values <- records$absolute_SMD[index[good]]
        q <- if (length(values)) as.numeric(quantile(values, c(0.25, 0.5, 0.75), type = 7, na.rm = FALSE)) else rep(NA_real_, 3L)
        data.frame(design = design, comparison_order = k, comparison_id = pair$comparison_id,
          group_A = pair$group_A, group_B = pair$group_B, covariate_order = j, covariate = covariate,
          median_absolute_SMD = q[2L], q25_absolute_SMD = q[1L], q75_absolute_SMD = q[3L],
          n_expected = length(ids), n_recorded = sum(present), n_valid = sum(good),
          n_unavailable = sum(!good), n_missing_records = sum(!present),
          n_flagged_invalid = sum(present & !good), quantile_type = 7L,
          summary_population = "valid_allocations_with_explicit_unavailable_denominator", row.names = NULL)
      }))
    }))
  }))
  rownames(result) <- NULL
  result
}
