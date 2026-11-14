test_that("intensity includes recurrent entries and paths starting down", {
  p <- 0.2
  q <- 0.3
  a <- c(0.4, 0.6)
  x <- MarkovChain(a, matrix(c(1 - p, p, q, 1 - q), 2, byrow = TRUE))
  time <- c(3, 1, 2, 1)
  expected <- p * (q / (p + q) + (a[1] - q / (p + q)) *
                    (1 - p - q)^(time - 1))
  expect_equal(dti(x, time, 1), expected)
  expect_equal(dti(x, 1, 1), ftp(x, 1, 1))
  expect_equal(vcov_dti(x, 1, 1), a[1]^2 * p * (1 - p) / 0.6)
  expect_gt(dti(x, 3, 1), ftp(x, 3, 1))
  down <- MarkovChain(c(0, 1), x$P)
  expect_equal(dti(down, 1:2, 1), c(0, q * p))
  expect_identical(vcov_dti(down, 1, 1), 0)
  # The q derivative is nonzero: down-row estimation error must contribute.
  expected_variance <- q^2 * p * (1 - p) / 0.6 +
    p^2 * q * (1 - q) / 0.4
  expect_equal(vcov_dti(down, 2, 1), expected_variance)
})

test_that("intensity variance matches independent stochastic derivatives", {
  P <- rbind(c(0.4, 0.1, 0.3, 0.2), c(0.2, 0.4, 0.1, 0.3),
             c(0.1, 0.2, 0.5, 0.2), c(0.3, 0.1, 0.2, 0.4))
  labels <- c("up-a", "down-a", "up-b", "down-b")
  dimnames(P) <- list(labels, labels)
  a <- setNames(c(0.4, 0.2, 0.3, 0.1), labels)
  x <- MarkovChain(a, P)
  balance <- t(P) - diag(4)
  balance[4, ] <- 1
  stationary <- solve(balance, c(0, 0, 0, 1))
  for (time in c(1, 2, 4)) {
    value <- function(B) {
      distribution <- a
      for (k in seq_len(time - 1)) distribution <- drop(distribution %*% B)
      sum(distribution[c(1, 3)] * rowSums(B[c(1, 3), c(2, 4)]))
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
    expect_equal(dti(x, time, c(1, 3)), value(P))
    expect_equal(vcov_dti(x, time, c(1, 3)), expected, tolerance = 1e-8)
    expect_equal(vcov_dti(x, time, labels[c(3, 1, 3)]), expected,
                 tolerance = 1e-8)
    permutation <- c(4, 3, 2, 1)
    permuted <- MarkovChain(a[permutation], P[permutation, permutation])
    expect_equal(vcov_dti(permuted, time, labels[c(1, 3)]), expected,
                 tolerance = 1e-8)
  }
})

test_that("intensity handles tiny exits, periodic and defective chains", {
  p <- 1e-18
  x <- MarkovChain(c(1, 0), matrix(c(1, p, 0.3, 0.7), 2, byrow = TRUE))
  expect_equal(dti(x, 1, 1) / p, 1)
  expect_equal(vcov_dti(x, 1, 1) / p, 1, tolerance = 1e-12)
  cycle <- MarkovChain(c(1, 0), matrix(c(0, 1, 1, 0), 2))
  expect_equal(dti(cycle, 1:4, 1), c(1, 0, 1, 0))
  expect_identical(vcov_dti(cycle, 4, 1), 0)
  P <- matrix(1 / 3, 3, 3) + 0.1 * outer(c(1, -1, 0), c(1, 1, -2))
  x <- MarkovChain(c(0.7, 0.2, 0.1), P)
  expect_equal(dti(x, 3:4, 1), rep(sum(P[1, 2:3]) / 3, 2))
  expect_true(is.finite(vcov_dti(x, 3, 1)))
})

test_that("intensity validates inputs and handles deterministic sets", {
  x <- MarkovChain(c(0.4, 0.6), matrix(c(0.8, 0.2, 0.3, 0.7), 2,
                                    byrow = TRUE))
  for (up in list(integer(), 1:2)) {
    expect_equal(dti(x, 1:3, up), rep(0, 3))
    expect_identical(vcov_dti(x, 2, up), 0)
  }
  reducible <- MarkovChain(c(1, 0), matrix(c(0.8, 0.2, 0, 1), 2, byrow = TRUE))
  expect_equal(dti(reducible, 1:3, 1), 0.2 * 0.8^(0:2))
  expect_error(vcov_dti(reducible, 1, 1), "irreducible")
  expect_error(vcov_dti(reducible, 1, integer()), "irreducible")
  single <- MarkovChain(1, matrix(1, 1, 1))
  expect_identical(dti(single, 1, 1), 0)
  expect_identical(vcov_dti(single, 1, 1), 0)
  for (f in list(dti, vcov_dti)) {
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
  expect_error(vcov_dti(x, 1:2, 1), "single")
})
