test_that("RG matches geometric rates and Corollary 2 at long horizons", {
  x <- MarkovChain(c(0.4, 0.6), matrix(c(0.8, 0.2, 0.3, 0.7), 2,
                                    byrow = TRUE))
  time <- c(3, 0, 1, 2, 1, 4000)
  expect_equal(rg(x, time, 1), -log(c(0.8, 0.4, rep(0.8, 4))))
  for (t in c(1, 2, 10, 4000)) {
    expect_equal(vcov_rg(x, t, 1), 0.2 / (0.8 * 0.6))
  }
  expect_identical(vcov_rg(x, 0, 1), 0)
})

test_that("RG retains tiny failure and survival probabilities", {
  p <- 1e-18
  rare <- MarkovChain(c(1, 0), matrix(c(1, p, 0.3, 0.7), 2, byrow = TRUE))
  expect_equal(rg(rare, 10, 1) / p, 1)
  expect_equal(vcov_rg(rare, 10, 1) / p, 1, tolerance = 1e-12)
  for (stay in c(1e-18, 1e-200)) {
    x <- MarkovChain(c(1, 0), matrix(c(stay, 1, 0.3, 0.7), 2, byrow = TRUE))
    expect_equal(rg(x, 1, 1), -log(stay))
    expected <- 1 / (stay * (0.3 / 1.3))
    expect_equal(vcov_rg(x, 1, 1) / expected, 1, tolerance = 1e-12)
  }
  x <- MarkovChain(c(1e-200, 1), matrix(c(0.8, 0.2, 0.3, 0.7), 2,
                                      byrow = TRUE))
  expect_equal(rg(x, 0, 1), -log(1e-200))
  expect_identical(vcov_rg(x, 0, 1), 0)
})

test_that("RG variance agrees with independent log-ratio derivatives", {
  P <- rbind(c(0.6, 0.3, 0.1), c(0.2, 0.5, 0.3), c(0.4, 0.1, 0.5))
  labels <- c("working", "failed", "degraded")
  dimnames(P) <- list(labels, labels)
  a <- setNames(c(0.7, 0.2, 0.1), labels)
  x <- MarkovChain(a, P)
  Q <- P[c(1, 3), c(1, 3)]
  for (time in c(1, 2, 5)) {
    value <- function(B) {
      power <- diag(2)
      for (k in seq_len(time - 1)) power <- power %*% B
      log(sum(a[c(1, 3)] %*% power) / sum(a[c(1, 3)] %*% power %*% B))
    }
    gradient <- matrix(0, 2, 2)
    h <- 1e-6
    for (i in seq_len(2)) {
      for (j in seq_len(2)) {
        H <- matrix(0, 2, 2)
        H[i, j] <- 1
        gradient[i, j] <- (value(Q + h * H) - value(Q - h * H)) / (2 * h)
      }
    }
    g <- as.vector(t(gradient))
    coordinates <- c(1, 3, 7, 9)
    expected <- drop(g %*% vcov_p(x)[coordinates, coordinates] %*% g)
    expect_equal(rg(x, time, c(1, 3)), value(Q))
    expect_equal(vcov_rg(x, time, c(1, 3)), expected, tolerance = 1e-8)
    expect_equal(vcov_rg(x, time, labels[c(3, 1, 3)]), expected,
                 tolerance = 1e-8)
  }
})

test_that("RG handles infinite and undefined rates explicitly", {
  x <- MarkovChain(c(1, 0), matrix(c(0, 1, 1, 0), 2))
  expect_equal(rg(x, 0:1, 1), c(0, Inf))
  expect_true(all(is.nan(rg(x, c(2, 5), 1))))
  expect_error(vcov_rg(x, 1, 1), "positive survival")
  expect_error(vcov_rg(x, 2, 1), "positive survival")
  expect_identical(rg(x, 0, integer()), Inf)
  expect_true(is.nan(rg(x, 1, integer())))
  expect_error(vcov_rg(x, 0, integer()), "positive survival")
  expect_equal(rg(x, 0:2, 1:2), rep(0, 3))
  expect_identical(vcov_rg(x, 2, 1:2), 0)
  single <- MarkovChain(1, matrix(1, 1, 1))
  expect_equal(rg(single, 0:2, 1), rep(0, 3))
  expect_identical(vcov_rg(single, 2, 1), 0)
})

test_that("RG validates inputs and variance irreducibility", {
  x <- MarkovChain(c(0.4, 0.6), matrix(c(0.8, 0.2, 0.3, 0.7), 2,
                                    byrow = TRUE))
  for (f in list(rg, vcov_rg)) {
    for (time in list(numeric(), -1, 1.5, NA, NaN, Inf, TRUE, "1",
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
  expect_error(vcov_rg(x, 1:2, 1), "single")
  reducible <- MarkovChain(c(1, 0), diag(2))
  expect_equal(rg(reducible, 3, 1), 0)
  expect_error(vcov_rg(reducible, 1, 1), "irreducible")
})
