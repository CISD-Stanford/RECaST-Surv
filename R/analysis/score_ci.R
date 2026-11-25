# A05 v001: deterministic approved score inversion; no inherited engine sourced.
valid_count <- function(x) is.numeric(x) && length(x) == 1L &&
  is.finite(x) && x >= 0 && x == floor(x)
valid_expected <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
pivot <- function(h, O, E) (O - h * E) / sqrt(h * E)

# Stable positive root of y^2 + z*y - O = 0, with h=y^2/E.
inverse_score <- function(z, O, E) {
  d <- sqrt(z * z + 4 * O)
  y <- ifelse(z >= 0, 2 * O / (d + z), (d - z) / 2)
  y * y / E
}

invert_cutoffs <- function(O, E, qL, qU) {
  out <- list(status = "unavailable", reason = "", lower = NA_real_,
    upper = NA_real_, length = NA_real_, zero_boundary = FALSE)
  reasons <- character()
  if (!valid_count(O)) reasons <- c(reasons, "invalid_final_O")
  if (!valid_expected(E)) reasons <- c(reasons, "nonpositive_or_nonfinite_final_E")
  if (length(qL) != 1L || length(qU) != 1L ||
      !all(is.finite(c(qL, qU))) || qL > qU)
    reasons <- c(reasons, "invalid_quantiles")
  if (length(reasons)) { out$reason <- paste(reasons, collapse = ";"); return(out) }
  if (O > 0) {
    limits <- c(inverse_score(qU, O, E), inverse_score(qL, O, E))
  } else if (qL >= 0) {
    out$status <- "empty"
    out$reason <- "zero_events_nonnegative_lower_cutoff"
    return(out)
  } else if (qU >= 0) {
    limits <- c(0, qL * qL / E)
    out$zero_boundary <- TRUE
  } else {
    limits <- c(qU * qU / E, qL * qL / E)
  }
  if (!all(is.finite(limits)) || any(limits < 0) || limits[1] > limits[2] ||
      (O > 0 && any(limits <= 0))) {
    out$reason <- "nonfinite_or_invalid_computed_endpoints"
    return(out)
  }
  out$status <- "nonempty"
  out$lower <- limits[1]; out$upper <- limits[2]
  out$length <- limits[2] - limits[1]
  out
}

construct_ci <- function(O, E, split_O, split_E) {
  q <- c(NA_real_, NA_real_)
  reasons <- character()
  if (length(split_O) != 100L || length(split_E) != 100L)
    reasons <- c(reasons, "requires_exactly_100_splits")
  if (!is.numeric(split_O) || any(!is.finite(split_O)) ||
      any(split_O < 0) || any(split_O != floor(split_O)))
    reasons <- c(reasons, "invalid_split_O")
  if (!is.numeric(split_E) || any(!is.finite(split_E)) || any(split_E <= 0))
    reasons <- c(reasons, "nonpositive_or_nonfinite_split_E")
  T <- rep(NA_real_, length(split_O))
  if (!length(reasons)) {
    T <- (split_O - split_E) / sqrt(split_E)
    if (!all(is.finite(T))) reasons <- c(reasons, "nonfinite_split_pivot")
    else q <- unname(quantile(T, c(.025, .975), type = 7, na.rm = FALSE))
  }
  out <- invert_cutoffs(O, E, q[1], q[2])
  if (length(reasons)) {
    # Invalid calibration never yields an interval, including for O=0.
    out$status <- "unavailable"
    out$reason <- paste(unique(c(reasons, strsplit(out$reason, ";", fixed = TRUE)[[1]])),
                        collapse = ";")
  }
  c(out, list(qL = q[1], qU = q[2], q_sum = sum(q),
              finite_split_pivots = sum(is.finite(T)), T = list(T)))
}

evaluate_reference <- function(ci, h_ref) {
  stopifnot(length(h_ref) == 1L, is.finite(h_ref), h_ref > 0)
  if (ci$status == "unavailable") return(c(contains = NA, below_lower = NA,
    above_upper = NA, empty_noncoverage = NA))
  if (ci$status == "empty") return(c(contains = FALSE, below_lower = FALSE,
    above_upper = FALSE, empty_noncoverage = TRUE))
  c(contains = ci$lower <= h_ref && h_ref <= ci$upper,
    below_lower = h_ref < ci$lower, above_upper = h_ref > ci$upper,
    empty_noncoverage = FALSE)
}

extract_method_splits <- function(splits, method) {
  stopifnot(method %in% c("median", "mean"))
  list(O = vapply(splits, function(s) s$methods[[method]]$O, numeric(1)),
       E = vapply(splits, function(s) s$methods[[method]]$E, numeric(1)),
       Z = vapply(splits, function(s) s$methods[[method]]$Z, numeric(1)))
}
