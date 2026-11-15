test_that("one-step intervals follow row-wise order and the multinomial variance", {
  states <- c("working", "failed")
  P <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE,
                      dimnames = list(states, states)))
  x <- MarkovChain(setNames(c(0.4, 0.6), states), P)
  result <- ci_p(x, n = 400, level = 0.9)
  expected <- c(0.8, 0.2, 0.3, 0.7)
  std_error <- sqrt(expected * (1 - expected) / rep(c(0.6, 0.4), each = 2) / 400)

  expect_s3_class(result, "data.frame")
  expect_named(result, c("from", "to", "estimate", "std.error", "lower", "upper"))
  expect_identical(result$from, rep(states, each = 2))
  expect_identical(result$to, rep(states, times = 2))
  expect_equal(result$estimate, expected)
  expect_equal(result$std.error, std_error)
  expect_equal(result$lower, expected - qnorm(0.95) * std_error)
  expect_equal(result$upper, expected + qnorm(0.95) * std_error)
  expect_identical(ci_p(P, n = 400, level = 0.9), result)
  expect_identical(ci_pk(x, k = 1, n = 400, level = 0.9), result)
})

test_that("three-step standard errors agree with independent finite differences", {
  P <- matrix(c(0.6, 0.3, 0.1, 0.15, 0.7, 0.15, 0.2, 0.3, 0.5),
              3, byrow = TRUE)
  stationary <- qr.solve(rbind(t(P) - diag(3), rep(1, 3)), c(0, 0, 0, 1))
  cube <- function(A) as.vector(t(A %*% A %*% A))
  variance <- numeric(9L)
  step <- 1e-5

  for (i in seq_len(nrow(P))) {
    derivative <- vapply(1:2, function(j) {
      H <- matrix(0, 3, 3)
      H[i, c(j, 3)] <- c(step, -step)
      (cube(P + H) - cube(P - H)) / (2 * step)
    }, numeric(9))
    p <- P[i, 1:2]
    row_covariance <- (diag(p) - tcrossprod(p)) / stationary[i]
    variance <- variance + rowSums((derivative %*% row_covariance) * derivative)
  }

  result <- ci_pk(StochMat(P), k = 3, n = 250)
  expect_equal(result$estimate, cube(P))
  expect_equal(result$std.error, sqrt(variance / 250), tolerance = 1e-8)
  expect_equal(result$lower, cube(P) - qnorm(0.975) * sqrt(variance / 250),
               tolerance = 1e-8)
  expect_equal(result$upper, cube(P) + qnorm(0.975) * sqrt(variance / 250),
               tolerance = 1e-8)
})

test_that("transition intervals handle deterministic entries and reject invalid inputs", {
  P <- StochMat(matrix(c(0, 1, 1, 0), 2, byrow = TRUE))
  result <- ci_pk(P, k = 2, n = 50)
  expect_equal(result$estimate, c(1, 0, 0, 1))
  expect_equal(result$std.error, rep(0, 4))
  expect_identical(result$lower, result$upper)

  singleton <- StochMat(matrix(1, 1, 1))
  expect_equal(ci_p(singleton, n = 10)$estimate, 1)
  expect_equal(ci_pk(singleton, k = 3, n = 10)$std.error, 0)

  for (k in list(0, -1, 1.5, NA, Inf, c(1, 2))) {
    expect_error(ci_pk(P, k = k, n = 50), "single positive integer")
  }
  expect_error(ci_p(P, n = 0), "n")
  expect_error(ci_pk(P, k = 2, n = 50, level = 1), "level")
  expect_error(ci_p(diag(2), n = 50), "MarkovChain.*StochMat")
  expect_error(ci_pk(StochMat(diag(2)), k = 2, n = 50), "irreducible")
})

test_that("deterministic k-step entries have zero standard error in a periodic chain", {
  P <- StochMat(rbind(
    c(0, 0.306620619678359, 0.644820231038023, 0.0485591492836182),
    c(1, 0, 0, 0),
    c(1, 0, 0, 0),
    c(1, 0, 0, 0)
  ))
  result <- ci_pk(P, k = 8, n = 50)

  expect_equal(result$estimate[1:4], c(1, 0, 0, 0))
  expect_equal(result$std.error[1:4], rep(0, 4))
  expect_true(all(result$std.error[c(6:8, 10:12, 14:16)] > 0))
})
