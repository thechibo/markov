test_that("failure-time probabilities and variances match two-state formulas", {
  p <- 0.2
  q <- 0.3
  x <- MarkovChain(c(0.4, 0.6), matrix(c(1 - p, p, q, 1 - q), 2,
                                    byrow = TRUE))
  time <- c(3, 1, 2, 1)
  expect_equal(ftp(x, time, 1), 0.4 * (1 - p)^(time - 1) * p)
  expect_equal(ftp(x, time, 1), rel(x, time - 1, 1) - rel(x, time, 1))
  for (t in c(1, 2, 5, 20)) {
    derivative <- if (t == 1) 0.4 else {
      0.4 * (1 - p)^(t - 2) * (1 - t * p)
    }
    expected <- derivative^2 * p * (1 - p) / (q / (p + q))
    expect_equal(vcov_ftp(x, t, 1), expected, tolerance = 1e-12)
  }
  # At p = 1/2 and t = 2 the two derivative terms cancel on stochastic rows.
  x <- MarkovChain(c(1, 0), matrix(c(0.5, 0.5, 0.3, 0.7), 2, byrow = TRUE))
  expect_identical(vcov_ftp(x, 2, 1), 0)
})

test_that("failure-time variance matches independent constrained derivatives", {
  P <- rbind(c(0.4, 0.1, 0.3, 0.2), c(0.2, 0.4, 0.1, 0.3),
             c(0.1, 0.2, 0.5, 0.2), c(0.3, 0.1, 0.2, 0.4))
  labels <- c("up-a", "down-a", "up-b", "down-b")
  dimnames(P) <- list(labels, labels)
  a <- setNames(c(0.4, 0.2, 0.3, 0.1), labels)
  x <- MarkovChain(a, P)
  U <- c(1, 3)
  balance <- t(P) - diag(4)
  balance[4, ] <- 1
  stationary <- solve(balance, c(0, 0, 0, 1))
  for (time in c(1, 2, 4)) {
    value <- function(B) {
      power <- diag(2)
      for (k in seq_len(time - 1)) power <- power %*% B[U, U]
      sum(a[U] %*% power %*% B[U, c(2, 4)])
    }
    gradient <- matrix(0, 4, 3)
    h <- 1e-6
    for (i in seq_len(4)) {
      for (j in seq_len(3)) {
        H <- matrix(0, 4, 4)
        H[i, j] <- 1
        H[i, 4] <- -1
        gradient[i, j] <- (value(P + h * H) - value(P - h * H)) / (2 * h)
      }
    }
    expected <- sum(vapply(seq_len(4), function(i) {
      p <- P[i, 1:3]
      drop(gradient[i, ] %*% (diag(p) - tcrossprod(p)) %*% gradient[i, ]) /
        stationary[i]
    }, numeric(1)))
    expect_equal(ftp(x, time, U), value(P))
    expect_equal(vcov_ftp(x, time, U), expected, tolerance = 1e-8)
    expect_equal(vcov_ftp(x, time, labels[c(3, 1, 3)]), expected,
                 tolerance = 1e-8)
  }
})

test_that("tiny exits and defective blocks do not need subtraction or eigenvectors", {
  p <- 1e-18
  q <- 0.3
  x <- MarkovChain(c(1, 0), matrix(c(1, p, q, 1 - q), 2, byrow = TRUE))
  expect_equal(ftp(x, 1, 1) / p, 1)
  expect_equal(vcov_ftp(x, 1, 1) / p, 1, tolerance = 1e-12)
  P <- rbind(c(0.5, 0.4, 0.1), c(0, 0.5, 0.5), c(0.2, 0.3, 0.5))
  x <- MarkovChain(c(1, 0, 0), P)
  expect_equal(ftp(x, 1:3, 1:2), c(0.1, 0.25, 0.225))
  expect_true(is.finite(vcov_ftp(x, 3, 1:2)))
})

test_that("failure-time APIs preserve boundary cases and reject malformed inputs", {
  x <- MarkovChain(c(0.4, 0.6), matrix(c(0.8, 0.2, 0.3, 0.7), 2,
                                    byrow = TRUE))
  for (up in list(integer(), 1:2)) {
    expect_equal(ftp(x, 1:3, up), rep(0, 3))
    expect_identical(vcov_ftp(x, 2, up), 0)
  }
  reducible <- MarkovChain(c(1, 0), matrix(c(0.8, 0.2, 0, 1), 2, byrow = TRUE))
  expect_equal(ftp(reducible, 1:3, 1), 0.2 * 0.8^(0:2))
  expect_error(vcov_ftp(reducible, 1, 1), "irreducible")
  expect_error(vcov_ftp(reducible, 1, integer()), "irreducible")
  single <- MarkovChain(1, matrix(1, 1, 1))
  expect_identical(ftp(single, 1, 1), 0)
  expect_identical(vcov_ftp(single, 1, 1), 0)
  down <- MarkovChain(c(0, 1), x$P)
  expect_identical(ftp(down, 2, 1), 0)
  expect_identical(vcov_ftp(down, 2, 1), 0)
  for (f in list(ftp, vcov_ftp)) {
    for (time in list(numeric(), 0, -1, 1.5, NA, NaN, Inf, TRUE, "1",
                     matrix(1), .Machine$integer.max)) {
      expect_error(f(x, time, 1), "time")
    }
    for (up in list(NA, Inf, 0, 3, 1.5, "unknown", matrix(1), TRUE)) {
      expect_error(f(x, 1, up), "up")
    }
    invalid <- x
    invalid$a[] <- c(2, 0)
    expect_error(f(invalid, 1, 1), "valid StochVec")
  }
  expect_error(vcov_ftp(x, 1:2, 1), "single")
})
