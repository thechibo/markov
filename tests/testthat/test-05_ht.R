test_that("test_prob returns a standard htest with correct transition counts", {
  x <- rep(c("a", "a", "b", "a", "b", "b"), 20)
  result <- test_prob(x, c("a", "b"), 0.5)
  expect_s3_class(result, "htest")
  expect_equal(result$observed, c(to_j = 40, other = 20))
  expect_equal(result$expected, c(to_j = 30, other = 30))
  expect_equal(result$estimate, c("p[a, b]" = 2 / 3))
  expect_equal(result$null.value, c("p[a, b]" = 0.5))
  expect_equal(result$statistic, c("X-squared" = 20 / 3))
  expect_equal(result$parameter, c(df = 1))
  expect_equal(result$p.value, pchisq(20 / 3, 1, lower.tail = FALSE))
  expect_identical(result$alternative, "two.sided")
  expect_identical(result$data.name, "x")
  expect_output(print(result), "Anderson-Goodman")
  reference <- prop.test(40, 60, p = 0.5, correct = FALSE)
  expect_equal(result$statistic, reference$statistic)
  expect_equal(result$p.value, reference$p.value)
})

test_that("test_prob aggregates all other destinations and excludes final origins", {
  x <- rep(c("a", "b", "a", "c", "a", "a"), 20)
  result <- test_prob(x, c("a", "b"), 0.3)
  expect_equal(result$observed, c(to_j = 20, other = 59))
  reference <- chisq.test(c(20, 59), p = c(0.3, 0.7))
  expect_equal(result$statistic, reference$statistic)
  self <- test_prob(x, c("a", "a"), 0.5)
  expect_equal(self$observed, c(to_j = 39, other = 40))
  expect_error(test_prob(c("a", "a", "b"), c("b", "a"), 0.5),
               "no observed departures")
})

test_that("test_prob matches state values and supports factors and numeric states", {
  x <- rep(c(10, 10, 20, 10, 20, 20), 20)
  result <- test_prob(x, c(10, 20), 0.5)
  expect_equal(unname(result$estimate), 2 / 3)
  expect_error(test_prob(x, c(1, 2), 0.5), "states observed")
  expect_equal(test_prob(factor(x), factor(c(10, 20)), 0.5)$statistic,
               result$statistic)
  close <- c(1, 1 + 1e-15)
  y <- rep(close[c(1, 1, 2, 1, 2, 2)], 20)
  expect_equal(test_prob(y, close, 0.5)$statistic, result$statistic)
  expect_equal(test_prob(rep("a", 21), c("a", "a"), 0.5)$statistic,
               c("X-squared" = 20))
})

test_that("test_prob validates trajectories, pairs and null probabilities", {
  x <- rep(c("a", "b"), 20)
  for (bad in list(NULL, 1)) {
    expect_error(test_prob(bad, c("a", "b"), 0.5), "length >= 2")
  }
  for (bad in list(list("a", "b"), array(x, c(2, 2, 10)))) {
    expect_error(test_prob(bad, c("a", "b"), 0.5), "atomic vector")
  }
  expect_error(test_prob(c("a", NA), c("a", "b"), 0.5), "NA values")
  for (bad in list(NULL, "a", c("a", "b", "a"), c("a", NA),
                   list("a", "b"), matrix(c("a", "b")))) {
    expect_error(test_prob(x, bad, 0.5), "exactly two")
  }
  expect_error(test_prob(x, c("a", "c"), 0.5), "states observed")
  for (bad in list(NULL, NA_real_, NaN, Inf, -0.1, 1.1, "0.5", TRUE,
                   c(0.2, 0.3), matrix(0.5), 0.5 + 0i)) {
    expect_error(test_prob(x, c("a", "b"), bad), "probability")
  }
  for (p0 in c(0, 1)) {
    expect_error(test_prob(x, c("a", "b"), p0), "0 < `p0` < 1", fixed = TRUE)
  }
  expect_warning(result <- test_prob(c("a", "b", "a"), c("a", "b"), 0.5),
                 "expected counts below 5")
  expect_true(is.finite(result$p.value))
})

test_that("test_prob pools only within-path transitions", {
  paths <- rbind(c("a", "a", "b"), c("b", "a", "b"))
  x <- paths[rep(1:2, each = 20), ]
  result <- test_prob(x, c("a", "b"), 0.5)
  expect_equal(result$observed, c(to_j = 40, other = 20))
  expect_equal(result$statistic, prop.test(40, 60, p = 0.5, correct = FALSE)$statistic)
  y <- rep(c("a", "a", "b", "a", "b", "b"), 20)
  expect_equal(test_prob(matrix(y, 1), c("a", "b"), 0.5)$statistic,
               test_prob(y, c("a", "b"), 0.5)$statistic)
  expect_error(test_prob(matrix("a", 2, 1), c("a", "a"), 0.5), "two columns")
})

test_that("test_const agrees with independent contingency-table calculations", {
  paths <- as.matrix(expand.grid(c("a", "b"), c("a", "b"), c("a", "b")))
  x <- paths[rep(1:8, times = c(12, 24, 36, 48, 60, 72, 84, 96)), ]
  result <- test_const(x)
  references <- lapply(c("a", "b"), function(state) {
    counts <- rbind(table(factor(x[x[, 1] == state, 2], levels = c("a", "b"))),
                    table(factor(x[x[, 2] == state, 3], levels = c("a", "b"))))
    chisq.test(counts, correct = FALSE)
  })
  expected_stat <- sum(vapply(references, function(z) unname(z$statistic), numeric(1)))
  expect_s3_class(result, "htest")
  expect_equal(unname(result$statistic), expected_stat)
  expect_equal(result$parameter, c(df = 2))
  expect_equal(result$p.value, pchisq(expected_stat, 2, lower.tail = FALSE))
  expect_equal(sum(result$observed), nrow(x) * 2)
  expect_equal(unname(rowSums(result$estimate)), c(1, 1))
  expect_equal(unname(result$estimate), unname(unclass(fit(x)$P)))
  expect_equal(result$expected[1, , ], t(unname(references[[1]]$expected)),
               ignore_attr = TRUE)
  expect_identical(result$data.name, "x")
  expect_output(print(result), "Anderson-Goodman")
  expect_equal(test_const(x[nrow(x):1, ])$statistic, result$statistic)
  balanced <- paths[rep(1:8, each = 20), ]
  expect_equal(test_const(balanced)$statistic, c("X-squared" = 0))
  expect_equal(test_const(balanced)$p.value, 1)
  # Four states and three transition times exercise non-binary dimensions.
  four <- as.matrix(expand.grid(1:4, 1:4, 1:4, 1:4))
  expect_equal(test_const(four)$parameter, c(df = 24))
  expect_equal(test_const(four)$statistic, c("X-squared" = 0))
})

test_that("test_const handles invalid inputs and inadequate counts explicitly", {
  for (bad in list(c("a", "b", "a"), matrix(1:3, 1), matrix(1:4, 2))) {
    expect_error(test_const(bad), "at least two paths and three observation times")
  }
  expect_error(test_const(matrix(1, 4, 3)), "two observed states")
  expect_error(test_const(rbind(c(1, 2, 1), c(1, 2, 1))), "positive origin counts")
  expect_error(test_const(rbind(c(1, 2, 1), c(2, 1, 2))), "positive pooled counts")
  expect_error(test_const(matrix(c(NA, rep(1, 5)), 2)), "NA values")
  expect_error(test_const(list(1, 2)), "atomic vector")
  paths <- as.matrix(expand.grid(1:2, 1:2, 1:2))
  expect_warning(result <- test_const(paths), "expected counts below 5")
  expect_equal(result$p.value, 1)
})
