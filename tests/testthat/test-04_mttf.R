test_that("MTTF and variance agree with a two-state analytic result", {
  p <- 0.2
  q <- 0.3
  x <- MarkovChain(c(1, 0), matrix(c(1 - p, p, q, 1 - q),
                                nrow = 2, byrow = TRUE))
  expected_variance <- (1 - p) * (p + q) / (q * p^3)
  expect_equal(mttf(x, 1), 1 / p)
  expect_equal(var_mttf(x, 1), expected_variance)
  expect_equal(var_mttf_sl(x, 1), expected_variance)
})

test_that("MTTF calculations preserve rare exit probabilities", {
  for (p in c(1e-12, 1e-16, 1e-18)) {
    q <- 0.5
    x <- MarkovChain(c(1, 0), matrix(c(1 - p, p, q, 1 - q),
                                  nrow = 2, byrow = TRUE))
    expected_variance <- (1 - p) * (p + q) / (q * p^3)
    expect_equal(mttf(x, 1) * p, 1, tolerance = 1e-12)
    expect_equal(var_mttf(x, 1) / expected_variance, 1, tolerance = 1e-12)
    expect_equal(var_mttf_sl(x, 1) / expected_variance, 1, tolerance = 1e-12)
  }
})

test_that("MTTF selectors use positions or names and ignore input order", {
  P <- matrix(c(0.6, 0.3, 0.1, 0.2, 0.5, 0.3, 0.4, 0.1, 0.5),
              nrow = 3, byrow = TRUE)
  x <- MarkovChain(c(0.7, 0.2, 0.1), P)
  # Numeric-looking names deliberately differ from their positions.
  labels <- c("3", "1", "2")
  dimnames(P) <- list(labels, labels)
  named <- MarkovChain(setNames(c(0.7, 0.2, 0.1), labels), P)
  for (f in list(mttf, var_mttf, var_mttf_sl, deriv_psi)) {
    for (up in list(c(1, 2), c(1, 3))) {
      expected <- f(x, up)
      expect_equal(f(x, rev(up)), expected)
      expect_equal(f(named, up), expected)
      expect_equal(f(named, labels[rev(up)]), expected)
    }
  }
})

test_that("MTTF variance agrees with independent constrained differences", {
  P <- matrix(c(0.6, 0.3, 0.1, 0.2, 0.5, 0.3, 0.4, 0.1, 0.5),
              nrow = 3, byrow = TRUE)
  a <- c(0.7, 0.2, 0.1)
  x <- MarkovChain(a, P)
  expect_equal(mttf(x, c(2, 1)), 34 / 7)
  gradient <- matrix(0, 3, 2)
  h <- 1e-6
  for (i in seq_len(3)) {
    for (j in seq_len(2)) {
      H <- matrix(0, 3, 3)
      H[i, j] <- 1
      H[i, 3] <- -1
      gradient[i, j] <- (
        mttf(MarkovChain(a, P + h * H), c(2, 1)) -
          mttf(MarkovChain(a, P - h * H), c(2, 1))
      ) / (2 * h)
    }
  }
  balance <- t(P) - diag(3)
  balance[3, ] <- 1
  stationary <- solve(balance, c(0, 0, 1))
  expected <- sum(vapply(seq_len(3), function(i) {
    covariance <- diag(P[i, 1:2]) - tcrossprod(P[i, 1:2])
    as.numeric(gradient[i, ] %*% covariance %*% gradient[i, ]) /
      stationary[i]
  }, numeric(1)))
  expect_equal(var_mttf(x, c(2, 1)), expected, tolerance = 1e-7)
  expect_equal(var_mttf_sl(x, c(2, 1)), expected, tolerance = 1e-7)
})

test_that("MTTF functions reject invalid state selections and objects", {
  x <- MarkovChain(c(0.7, 0.2, 0.1), matrix(
    c(0.6, 0.3, 0.1, 0.2, 0.5, 0.3, 0.4, 0.1, 0.5), 3, byrow = TRUE
  ))
  for (f in list(mttf, var_mttf, var_mttf_sl, deriv_psi)) {
    expect_error(f(x, c(1, 1)), "distinct")
    expect_error(f(x, c("1", "1")), "distinct")
    expect_error(f(x, integer()), "nonempty")
    expect_error(f(x, NA_integer_), "nonempty")
    expect_error(f(x, TRUE), "nonempty")
    expect_error(f(x, c(1, 4)), "indices")
    expect_error(f(x, 1.5), "indices")
    expect_error(f(x, "unknown"), "unknown")
    expect_error(f(x, 1:3), "proper subset")
    expect_error(f(unclass(x), 1), "MarkovChain")
    invalid <- x
    invalid$a[1] <- -1
    expect_error(f(invalid, 1), "valid StochVec")
    invalid <- x
    invalid$P[1, ] <- c(-0.1, 0.5, 0.6)
    expect_error(f(invalid, 1), "valid StochMat")
  }
})
