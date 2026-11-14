test_that("stationary availability matches the two-state analytic result", {
  for (rates in list(c(0.2, 0.3), c(1e-8, 2e-8))) {
    p <- rates[1]
    q <- rates[2]
    x <- MarkovChain(c(1, 0), matrix(c(1 - p, p, q, 1 - q), 2, byrow = TRUE))
    pi <- c(q, p) / (p + q)
    expected <- (q / (p + q)^2)^2 * p * (1 - p) / pi[1] +
      (p / (p + q)^2)^2 * q * (1 - q) / pi[2]
    expect_equal(savail(x, 1), pi[1])
    expect_equal(savail(x, 2), pi[2])
    expect_equal(vcov_savail(x, 1) / expected, 1, tolerance = 1e-7)
    expect_equal(vcov_savail(x, 2) / expected, 1, tolerance = 1e-7)
    other <- MarkovChain(c(0, 1), x$P)
    expect_equal(savail(other, 1), savail(x, 1))
    expect_equal(vcov_savail(other, 1), vcov_savail(x, 1))
  }
})

test_that("stationary availability covariance includes cross-state terms", {
  P <- rbind(c(0.4, 0.1, 0.3, 0.2), c(0.2, 0.4, 0.1, 0.3),
             c(0.1, 0.2, 0.5, 0.2), c(0.3, 0.1, 0.2, 0.4))
  labels <- c("up-a", "down-a", "up-b", "down-b")
  dimnames(P) <- list(labels, labels)
  a <- setNames(c(0.4, 0.2, 0.3, 0.1), labels)
  x <- MarkovChain(a, P)
  stationary <- function(B) {
    balance <- t(B) - diag(4)
    balance[4, ] <- 1
    solve(balance, c(0, 0, 0, 1))
  }
  pi <- stationary(P)
  for (up in list(c(1, 3), c(1, 2, 3))) {
    value <- function(B) sum(stationary(B)[up])
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
      drop(gradient[i, ] %*% (diag(p) - tcrossprod(p)) %*% gradient[i, ]) / pi[i]
    }, numeric(1)))
    expect_equal(savail(x, up), value(P))
    expect_equal(vcov_savail(x, up), expected, tolerance = 1e-8)
    expect_equal(vcov_savail(x, setdiff(1:4, up)), expected, tolerance = 1e-8)
    expect_equal(vcov_savail(x, labels[c(rev(up), up[1])]), expected,
                 tolerance = 1e-8)
    permutation <- c(4, 3, 2, 1)
    permuted <- MarkovChain(a[permutation], P[permutation, permutation])
    expect_equal(savail(permuted, labels[up]), value(P))
    expect_equal(vcov_savail(permuted, labels[up]), expected, tolerance = 1e-8)
    Sigma <- vcov_stat(x)
    expect_equal(vcov_savail(x, up), sum(Sigma[up, up, drop = FALSE]))
  }
})

test_that("stationary availability handles periodic and deterministic cases", {
  x <- MarkovChain(c(1, 0), matrix(c(0, 1, 1, 0), 2))
  expect_equal(savail(x, 1), 0.5)
  expect_identical(vcov_savail(x, 1), 0)
  expect_identical(savail(x, integer()), 0)
  expect_identical(savail(x, 1:2), 1)
  expect_identical(vcov_savail(x, integer()), 0)
  expect_identical(vcov_savail(x, 1:2), 0)
  expect_equal(savail(x, c(1, 1)), 0.5)
  single <- MarkovChain(1, matrix(1, 1, 1))
  expect_identical(savail(single, 1), 1)
  expect_identical(vcov_savail(single, 1), 0)
})

test_that("stationary availability enforces input and irreducibility contracts", {
  x <- MarkovChain(c(0.4, 0.6), matrix(c(0.8, 0.2, 0.3, 0.7), 2,
                                    byrow = TRUE))
  reducible <- MarkovChain(c(1, 0), matrix(c(0.8, 0.2, 0, 1), 2, byrow = TRUE))
  multiple <- MarkovChain(c(1, 0), diag(2))
  for (f in list(savail, vcov_savail)) {
    for (up in list(NA, Inf, 0, 3, 1.5, "unknown", matrix(1), TRUE)) {
      expect_error(f(x, up), "up")
    }
    for (up in list(1, integer(), 1:2)) {
      expect_error(f(reducible, up), "irreducible")
      expect_error(f(multiple, up), "irreducible")
    }
    invalid <- x
    invalid$a[] <- c(2, 0)
    expect_error(f(invalid, 1), "valid StochVec")
  }
})
