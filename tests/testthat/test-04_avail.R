test_that("availability agrees with the two-state formula", {
  p <- 0.2
  q <- 0.3
  a <- c(0.4, 0.6)
  x <- MarkovChain(a, matrix(c(1 - p, p, q, 1 - q), 2, byrow = TRUE))
  time <- c(3L, 1L, 2L, 1L)
  expected <- q / (p + q) + (a[1] - q / (p + q)) * (1 - p - q)^time
  expect_equal(avail(x, time, 1), expected)
  expect_equal(avail(x, time, c(1, 1)), expected)
  expect_equal(avail(x, time, 2), 1 - expected)
  expect_equal(avail(x, time, 1:2), rep(1, length(time)))
})

test_that("availability variance agrees with an independent delta method", {
  p <- 0.2
  q <- 0.3
  a <- c(0.4, 0.6)
  x <- MarkovChain(a, matrix(c(1 - p, p, q, 1 - q), 2, byrow = TRUE))
  stationary <- c(q, p) / (p + q)
  parameter_variance <- c(p * (1 - p), q * (1 - q)) / stationary
  availability <- function(p, q, time) {
    q / (p + q) + (a[1] - q / (p + q)) * (1 - p - q)^time
  }
  h <- 1e-6
  for (time in c(1L, 3L)) {
    gradient <- c(
      (availability(p + h, q, time) - availability(p - h, q, time)) / (2 * h),
      (availability(p, q + h, time) - availability(p, q - h, time)) / (2 * h)
    )
    expected <- sum(gradient^2 * parameter_variance)
    expect_equal(vcov_avail(x, time, 1), expected, tolerance = 1e-7)
    expect_equal(vcov_avail(x, time, 2), expected, tolerance = 1e-7)
    expect_equal(vcov_avail_sl(x, time, 1), expected, tolerance = 1e-7)
    expect_equal(vcov_avail_sl(x, time, 2), expected, tolerance = 1e-7)
  }
})

test_that("the availability reference includes every initial state", {
  a <- c(0.4, 0.6)
  x <- MarkovChain(a, matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE))
  expected <- a[1]^2 * 0.8 * 0.2 / 0.6 + a[2]^2 * 0.3 * 0.7 / 0.4
  expect_equal(vcov_avail_sl(x, 1, 1), expected)
  expect_equal(vcov_avail_sl(x, 1, 2), expected)
})

# test_that("rare transitions retain positive variance for complementary up sets", {
#   entry <- 1e-12
#   stay <- 1e-18
#   P <- rbind(c(1 - entry, entry, 0),
#              c(0.45, stay, 0.55 - stay),
#              c(0.55, 0, 0.45))
#   x <- MarkovChain(c(0, 1, 0), P)
#   # From stationary flow into states 2 and 3, relative to state 1.
#   weights <- c(1, entry / (1 - stay),
#                entry / (1 - stay) * (0.55 - stay) / 0.55)
#   expected <- stay * (1 - stay) * sum(weights) / weights[2]
#   for (f in list(vcov_avail, vcov_avail_sl)) {
#     expect_equal(f(x, 1L, c(1L, 3L)) / expected, 1, tolerance = 1e-12)
#     expect_equal(f(x, 1L, 2L) / expected, 1, tolerance = 1e-12)
#   }
# })

# test_that("availability variance handles deterministic up sets exactly", {
#   entry <- 1e-12
#   P <- rbind(c(1 - entry, entry, 0),
#              c(0.45, 0.2, 0.35),
#              c(0.55, 0.2, 0.25))
#   x <- MarkovChain(c(0, 1, 0), P)
#   reducible <- MarkovChain(c(0, 1, 0), diag(3))
#   for (f in list(vcov_avail, vcov_avail_sl)) {
#     expect_identical(f(x, 3L, 1:3), 0)
#     expect_identical(f(x, 3L, integer()), 0)
#     expect_error(f(reducible, 3L, 1:3), "irreducible")
#     expect_error(f(reducible, 3L, integer()), "irreducible")
#   }
# })

test_that("availability supports ordinary chains with complex eigenvalues", {
  P <- matrix(c(0.6, 0.3, 0.1, 0.2, 0.5, 0.3, 0.4, 0.1, 0.5),
              3, byrow = TRUE)
  a <- c(0.7, 0.2, 0.1)
  x <- MarkovChain(a, P)
  # Compute stationary probabilities independently of package helpers.
  balance <- t(P) - diag(3)
  balance[3, ] <- 1
  stationary <- solve(balance, c(0, 0, 1))
  expected <- sum(a^2 * P[, 1] * (1 - P[, 1]) / stationary)
  expect_equal(vcov_avail(x, 1, 1), expected)
  expect_no_warning(expect_equal(avail(x, 1:3, 1), c(0.5, 0.436, 0.4216)))
  expect_equal(vcov_avail(x, 3, 1), vcov_avail_sl(x, 3, 1))
})

test_that("availability does not require independent eigenvectors", {
  P <- matrix(c(0.5, 0.5, 0, 0, 0.5, 0.5, 0, 0, 1), 3, byrow = TRUE)
  x <- MarkovChain(c(1, 0, 0), P)
  expect_equal(avail(x, c(3, 1, 2), 2), c(0.375, 0.5, 0.5))

  # A positive defective chain: P^2 is the uniform transition matrix.
  P <- matrix(1 / 3, 3, 3) + 0.1 * outer(c(1, -1, 0), c(1, 1, -2))
  x <- MarkovChain(c(0.7, 0.2, 0.1), P)
  expect_equal(avail(x, 1:3, 1), c(1 / 3 + 0.05, 1 / 3, 1 / 3))
  expect_equal(vcov_avail(x, 3, 1), vcov_avail_sl(x, 3, 1))
})

test_that("availability resolves numeric selectors with named states", {
  labels <- c("working", "failed")
  P <- matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE,
              dimnames = list(labels, labels))
  x <- MarkovChain(setNames(c(0.4, 0.6), labels), P)
  expect_equal(avail(x, 1, 1), avail(x, 1, "working"))
  expect_equal(vcov_avail(x, 1, 1), vcov_avail(x, 1, "working"))
  expect_equal(vcov_avail_sl(x, 1, 1), vcov_avail_sl(x, 1, "working"))
})

test_that("availability validates its object and scalar variance time", {
  x <- MarkovChain(c(0.4, 0.6), matrix(c(0.8, 0.2, 0.3, 0.7),
                                    2, byrow = TRUE))
  expect_error(vcov_avail(x, c(1, 2), 1), "time")
  expect_error(avail(x, -1, 1), "time")
  invalid <- x
  invalid$a[] <- c(2, 0)
  expect_error(avail(invalid, 1, 1), "valid StochVec")
  for (f in list(vcov_avail, vcov_avail_sl)) {
    expect_error(f(x, c(1, 2), 1), "single positive integer")
    expect_error(f(invalid, 1, 1), "valid StochVec")
    expect_error(f(x, 1, 1.5), "indices")
    expect_error(f(x, 1, NA_integer_), "up")
  }
})

test_that("the component-wise reference handles noncontiguous named up states", {
  P <- matrix(c(0.6, 0.3, 0.1, 0.2, 0.5, 0.3, 0.4, 0.1, 0.5),
              3, byrow = TRUE)
  labels <- c("third", "first", "second")
  dimnames(P) <- list(labels, labels)
  a <- setNames(c(0.7, 0.2, 0.1), labels)
  x <- MarkovChain(a, P)
  # Independent stochastic perturbations, with column 2 compensating each row.
  value <- function(Q) sum((a %*% Q %*% Q %*% Q)[c(1, 3)])
  gradient <- matrix(0, 3, 2)
  h <- 1e-6
  for (i in seq_len(3)) {
    for (j in seq_len(2)) {
      H <- matrix(0, 3, 3)
      H[i, c(1, 3)[j]] <- 1
      H[i, 2] <- -1
      gradient[i, j] <- (value(P + h * H) - value(P - h * H)) / (2 * h)
    }
  }
  balance <- t(P) - diag(3)
  balance[3, ] <- 1
  stationary <- solve(balance, c(0, 0, 1))
  expected <- sum(vapply(seq_len(3), function(i) {
    p <- P[i, c(1, 3)]
    drop(gradient[i, ] %*% (diag(p) - tcrossprod(p)) %*% gradient[i, ]) /
      stationary[i]
  }, numeric(1)))
  expect_equal(vcov_avail_sl(x, 3, c(3, 1)), expected, tolerance = 1e-7)
  expect_equal(vcov_avail(x, 3, c(3, 1)), expected, tolerance = 1e-7)
  expect_equal(vcov_avail_sl(x, 3, labels[c(3, 1, 3)]), expected,
               tolerance = 1e-7)
})
