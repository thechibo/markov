ci_test_chain <- function() {
  MarkovChain(c(0.4, 0.6), matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE))
}

test_that("all reliability CIs use the corresponding point estimate and variance", {
  x <- ci_test_chain()
  variances <- c(rel = "vcov_rel", avail = "vcov_avail", dssc = "vcov_dssc",
                 main = "vcov_main", ftp = "vcov_ftp", dti = "vcov_dti",
                 bmpf = "vcov_bmp", rgf = "vcov_rg", savail = "vcov_savail",
                 mttf = "var_mttf", mttr = "var_mttr")
  for (statistic in names(variances)) {
    timed <- !statistic %in% c("savail", "mttf", "mttr")
    args <- list(x = x)
    if (timed) args$time <- c(3L, 1L, 3L)
    states <- if (statistic %in% c("dssc", "main", "mttr")) "down" else "up"
    args[[states]] <- 1L
    ci <- get(paste0("ci_", statistic))
    result <- do.call(ci, c(args, list(n = 400, level = 0.9)))
    estimate <- do.call(get(statistic), args)
    variance <- if (timed) vapply(args$time, function(t) {
      a <- args
      a$time <- t
      do.call(get(variances[[statistic]]), a)
    }, numeric(1)) else do.call(get(variances[[statistic]]), args)
    expect_s3_class(result, "data.frame")
    expect_named(result, c("estimate", "std.error", "lower", "upper"))
    expect_equal(result$estimate, estimate)
    expect_equal(result$std.error, sqrt(variance / 400))
    expect_equal(result$lower, estimate - qnorm(0.95) * sqrt(variance / 400))
    expect_equal(result$upper, estimate + qnorm(0.95) * sqrt(variance / 400))
    larger <- do.call(ci, c(args, list(n = 1600, level = 0.9)))
    expect_equal(larger$std.error, result$std.error / 2)
    wider <- do.call(ci, c(args, list(n = 400, level = 0.99)))
    expect_true(all(wider$lower < result$lower))
    expect_true(all(wider$upper > result$upper))
    scalar_args <- args
    if (timed) scalar_args$time <- 3L
    expect_equal(nrow(do.call(ci, c(scalar_args, list(n = 400)))), 1L)
  }
})

test_that("two-state CIs agree with independent analytic variances", {
  x <- ci_test_chain()
  # p = P11, q = P21, pi1 = q / (1 - p + q).
  p <- 0.8
  q <- 0.3
  a <- 0.4
  pi1 <- q / (1 - p + q)
  vp <- p * (1 - p) / pi1
  vq <- q * (1 - q) / (1 - pi1)
  time <- c(0L, 1L, 4L)
  expect_equal(ci_rel(x, time, 1, 100)$std.error,
               a * time * p^pmax(time - 1, 0) * sqrt(vp / 100))
  expect_equal(ci_mttf(x, 1, 100)$std.error,
               a / (1 - p)^2 * sqrt(vp / 100))
  expect_equal(ci_mttr(x, 2, 100)$std.error,
               (1 - a) / q^2 * sqrt(vq / 100))
  expect_equal(ci_bmpf(x, 1:3, 1, 100)$std.error, rep(sqrt(vp / 100), 3))
  expect_equal(ci_rgf(x, 1:3, 1, 100)$std.error, rep(sqrt(vp / 100) / p, 3))
  stationary_var <- (q^2 * vp + (1 - p)^2 * vq) / (1 - p + q)^4
  expect_equal(ci_savail(x, 1, 100)$std.error, sqrt(stationary_var / 100))
})

test_that("CI controls reject invalid values before computing statistics", {
  invalid_n <- list(NULL, numeric(), c(1, 2), NA_real_, NaN, Inf, -Inf,
                    0, -1, 1.5, "100", TRUE, matrix(100))
  invalid_level <- list(NULL, numeric(), c(0.9, 0.95), NA_real_, NaN, Inf,
                        -Inf, 0, 1, -0.1, 1.1, "0.95", TRUE, matrix(0.95))
  for (name in c("rel", "avail", "dssc", "main", "ftp", "dti", "bmpf",
                 "rgf", "savail", "mttf", "mttr")) {
    ci <- get(paste0("ci_", name))
    for (n in invalid_n) expect_error(ci(n = n), "`n`", fixed = TRUE)
    for (level in invalid_level) {
      expect_error(ci(n = 100, level = level), "`level`", fixed = TRUE)
    }
  }
})

test_that("degenerate intervals and inherited domain restrictions are explicit", {
  x <- ci_test_chain()
  for (name in c("rel", "dssc", "main", "bmpf", "rgf")) {
    result <- get(paste0("ci_", name))(x, 0L, 1, 100)
    expect_equal(result$std.error, 0)
    expect_equal(result$lower, result$estimate)
    expect_equal(result$upper, result$estimate)
  }
  expect_equal(ci_savail(x, integer(), 100)$lower, 0)
  expect_equal(ci_rel(x, 1:3, 1:2, 100)$upper, rep(1, 3))
  expect_equal(ci_main(x, 1:3, 2, 100)$lower,
               1 - ci_dssc(x, 1:3, 2, 100)$upper)
  expect_error(ci_rel(x, numeric(), 1, 100), "time")
  expect_error(ci_avail(x, 0, 1, 100), "time")
  expect_error(ci_ftp(x, 0, 1, 100), "time")
  expect_error(ci_rel(x, 1, "missing", 100), "unknown")
  expect_error(ci_mttf(x, 1:2, 100), "proper subset")
  expect_error(ci_bmpf(x, 1, integer(), 100), "survival")
  expect_error(ci_rgf(x, 0, integer(), 100), "survival")
  reducible <- MarkovChain(c(1, 0), diag(2))
  expect_error(ci_rel(reducible, 1, 1, 100), "irreducible")
  expect_lt(ci_rel(x, 1, 1, 1)$lower, 0) # No clipping.
})

test_that("Normal-approximation arithmetic handles extreme scales without premature underflow", {
  expect_error(norm_ci(Inf, 1, 100, 0.95), "finite point")
  for (v in c(-1, NA, NaN, Inf)) {
    expect_error(norm_ci(1, v, 100, 0.95), "finite nonnegative")
  }
  expect_equal(norm_ci(0, 1e-300, 1e100, 0.95)$std.error / 1e-200, 1)
  result <- ci_rel(ci_test_chain(), 1, 1, 100, level = 1 - .Machine$double.eps / 2)
  expect_true(all(is.finite(unlist(result))))
})
