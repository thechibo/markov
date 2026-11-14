test_that("reliability excludes repairs and retains initial down mass", {
  x <- MarkovChain(c(0.4, 0.6), matrix(c(0.8, 0.2, 0.3, 0.7), 2,
                                    byrow = TRUE))
  time <- c(3, 0, 1, 2, 1)
  expect_equal(rel(x, time, 1), 0.4 * 0.8^time)
  expect_lt(rel(x, 2, 1), avail(x, 2, 1))
  for (f in list(vcov_rel, vcov_rel_sl)) {
    for (t in c(1, 2, 5, 50)) {
      expected <- (0.4 * t * 0.8^(t - 1))^2 * 0.8 * 0.2 / 0.6
      expect_equal(f(x, t, 1), expected)
    }
  }
})

test_that("reliability variance matches independent stochastic perturbations", {
  P <- rbind(c(0.6, 0.3, 0.1), c(0.2, 0.5, 0.3), c(0.4, 0.1, 0.5))
  labels <- c("third", "first", "second")
  dimnames(P) <- list(labels, labels)
  a <- setNames(c(0.7, 0.2, 0.1), labels)
  x <- MarkovChain(a, P)
  U <- c(1, 3)
  balance <- t(P) - diag(3)
  balance[3, ] <- 1
  stationary <- solve(balance, c(0, 0, 1))
  for (time in c(1, 2, 4)) {
    value <- function(B) {
      Q <- B[U, U]
      power <- diag(2)
      for (k in seq_len(time)) power <- power %*% Q
      sum(a[U] %*% power)
    }
    gradient <- matrix(0, 2, 2)
    h <- 1e-6
    for (i in seq_len(2)) {
      for (j in seq_len(2)) {
        H <- matrix(0, 3, 3)
        H[U[i], U[j]] <- 1
        H[U[i], 2] <- -1
        gradient[i, j] <- (value(P + h * H) - value(P - h * H)) / (2 * h)
      }
    }
    expected <- sum(vapply(seq_len(2), function(i) {
      p <- P[U[i], U]
      drop(gradient[i, ] %*% (diag(p) - tcrossprod(p)) %*% gradient[i, ]) /
        stationary[U[i]]
    }, numeric(1)))
    expect_equal(rel(x, time, U), value(P))
    for (f in list(vcov_rel, vcov_rel_sl)) {
      expect_equal(f(x, time, U), expected, tolerance = 1e-8)
      expect_equal(f(x, time, labels[c(3, 1, 3)]), expected,
                   tolerance = 1e-8)
    }
  }
})

test_that("the Kronecker contraction matches the full Equation 52 operator", {
  P <- matrix(c(0.4, 0.1, 0.3, 0.2, 0.2, 0.4, 0.1, 0.3,
                0.1, 0.2, 0.5, 0.2, 0.3, 0.1, 0.2, 0.4), 4, byrow = TRUE)
  x <- MarkovChain(c(0.1, 0.2, 0.3, 0.4), P)
  U <- c(1, 3, 4)
  Q <- P[U, U]
  powers <- list(diag(3), Q, Q %*% Q)
  M <- Reduce(`+`, lapply(1:3, function(i) {
    kronecker(powers[[i]], t(powers[[4 - i]]))
  }))
  C <- kronecker(matrix(x$a[U], 1), matrix(1, 1, 3))
  coordinates <- c(1, 3, 4, 9, 11, 12, 13, 15, 16)
  Sigma <- vcov_p(x)[coordinates, coordinates]
  expected <- drop(C %*% M %*% Sigma %*% t(M) %*% t(C))
  expect_equal(vcov_rel(x, 3, U), expected)
  expect_equal(vcov_rel_sl(x, 3, U), expected)
})

test_that("reliability handles defective and nearly reducible chains", {
  P <- rbind(c(0.5, 0.4, 0.1), c(0, 0.5, 0.5), c(0.2, 0.3, 0.5))
  x <- MarkovChain(c(1, 0, 0), P)
  # The operational block is a Jordan block.
  expect_equal(rel(x, 1:5, 1:2), 0.5^(1:5) + (1:5) * 0.4 * 0.5^(0:4))
  expect_equal(vcov_rel(x, 5, 1:2), vcov_rel_sl(x, 5, 1:2))
  p <- 1e-10
  q <- 2e-10
  x <- MarkovChain(c(1, 0), matrix(c(1 - p, p, q, 1 - q), 2, byrow = TRUE))
  expected <- (10 * (1 - p)^9)^2 * (1 - p) * p / (q / (p + q))
  for (f in list(vcov_rel, vcov_rel_sl)) {
    expect_equal(f(x, 10, 1) / expected, 1, tolerance = 1e-7)
  }
})

test_that("deterministic reliability cases and irreducibility are explicit", {
  x <- MarkovChain(c(0.4, 0.6), matrix(c(0.8, 0.2, 0.3, 0.7), 2,
                                    byrow = TRUE))
  expect_equal(rel(x, 0:2, integer()), rep(0, 3))
  expect_equal(rel(x, 0:2, 1:2), rep(1, 3))
  expect_equal(rel(x, 0, c(1, 1)), 0.4)
  down <- MarkovChain(c(0, 1), x$P)
  reducible <- MarkovChain(c(1, 0), diag(2))
  expect_equal(rel(reducible, 0:2, 1), rep(1, 3))
  for (f in list(vcov_rel, vcov_rel_sl)) {
    expect_identical(f(x, 0, 1), 0)
    expect_identical(f(x, 3, integer()), 0)
    expect_identical(f(x, 3, 1:2), 0)
    expect_identical(f(down, 3, 1), 0)
    expect_error(f(reducible, 0, 1), "irreducible")
    expect_error(f(reducible, 3, integer()), "irreducible")
  }
  single <- MarkovChain(1, matrix(1, 1, 1))
  expect_equal(rel(single, 0:2, 1), rep(1, 3))
  expect_identical(vcov_rel(single, 2, 1), 0)
  expect_identical(vcov_rel_sl(single, 2, 1), 0)
})

test_that("reliability rejects malformed inputs", {
  x <- MarkovChain(c(0.4, 0.6), matrix(c(0.8, 0.2, 0.3, 0.7), 2,
                                    byrow = TRUE))
  for (f in list(rel, vcov_rel, vcov_rel_sl)) {
    for (time in list(numeric(), -1, 1.5, NA, NaN, Inf, "1", TRUE,
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
  expect_error(vcov_rel(x, 1:2, 1), "single")
  expect_error(vcov_rel_sl(x, 1:2, 1), "single")
})
