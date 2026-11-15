# Shared validation and arithmetic for pointwise normal-approximation intervals.
validate_ci_parameters <- function(n, level) {
  if (!is.numeric(n) || !is.null(dim(n)) || length(n) != 1L ||
      is.na(n) || !is.finite(n) || n <= 0 || n != floor(n)) {
    stop("`n` must be a single finite positive integer.", call. = FALSE)
  }
  if (!is.numeric(level) || !is.null(dim(level)) || length(level) != 1L ||
      is.na(level) || !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be a single finite number strictly between 0 and 1.",
         call. = FALSE)
  }
  invisible(NULL)
}

norm_ci <- function(estimate, variance, n, level) {
  if (any(!is.finite(estimate))) {
    stop("Confidence intervals require finite point estimates.", call. = FALSE)
  }
  if (any(!is.finite(variance)) || any(variance < 0)) {
    stop("Confidence intervals require finite nonnegative variances.",
         call. = FALSE)
  }
  # Separate square roots avoid underflow in variance / n. The upper-tail
  # quantile avoids rounding (1 + level) / 2 to one near level = 1.
  std.error <- sqrt(variance) / sqrt(n)
  margin <- stats::qnorm((1 - level) / 2, lower.tail = FALSE) * std.error
  data.frame(estimate = as.numeric(estimate), std.error = std.error,
             lower = estimate - margin, upper = estimate + margin,
             row.names = NULL)
}
