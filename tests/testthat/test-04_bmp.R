test_that("BMP matches the geometric hazard without underflow", {
  x <- MarkovChain(c(0.4, 0.6), matrix(c(0.8, 0.2, 0.3, 0.7), 2,
                                    byrow = TRUE))
  expect_equal(bmp(x, c(3, 0, 1, 2, 1, 4000), 1), c(0.2, 0.6, rep(0.2, 4)))
  for (t in c(1, 2, 10, 4000)) {
    expect_equal(vcov_bmp(x, t, 1), 0.8 * 0.2 / 0.6)
  }
  expect_identical(vcov_bmp(x, 0, 1), 0)
  p <- 1e-18
  rare <- MarkovChain(c(1, 0), matrix(c(1, p, 0.3, 0.7), 2, byrow = TRUE))
  expect_equal(bmp(rare, 10, 1) / p, 1)
  expect_equal(vcov_bmp(rare, 10, 1) / p, 1, tolerance = 1e-12)
})

test_that("BMP agrees with Theorem 7 quotient derivatives", {
  P <- rbind(c(0.4, 0.1, 0.3, 0.2), c(0.2, 0.4, 0.1, 0.3),
             c(0.1, 0.2, 0.5, 0.2), c(0.3, 0.1, 0.2, 0.4))
  labels <- c("up-a", "down-a", "up-b", "down-b")
  dimnames(P) <- list(labels, labels)
  a <- setNames(c(0.4, 0.2, 0.3, 0.1), labels)
  x <- MarkovChain(a, P)
  U <- c(1, 3)
  Q <- P[U, U]
  for (time in c(1, 2, 5)) {
    value <- function(B) {
      power <- diag(2)
      for (k in seq_len(time - 1)) power <- power %*% B
      1 - sum(a[U] %*% power %*% B) / sum(a[U] %*% power)
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
    # Restricted row-wise coordinates in the full transition covariance.
    g <- as.vector(t(gradient))
    coordinates <- c(1, 3, 9, 11)
    expected <- drop(g %*% vcov_p(x)[coordinates, coordinates] %*% g)
    expect_equal(bmp(x, time, U), value(Q))
    expect_equal(vcov_bmp(x, time, U), expected, tolerance = 1e-8)
    expect_equal(vcov_bmp(x, time, labels[c(3, 1, 3)]), expected,
                 tolerance = 1e-8)
    expect_equal(bmp(x, time, U), ftp(x, time, U) / rel(x, time - 1, U))
    # Scaling the initial operational mass leaves the conditional law fixed.
    conditional <- MarkovChain(setNames(c(4 / 7, 0, 3 / 7, 0), labels), P)
    expect_equal(vcov_bmp(conditional, time, U), vcov_bmp(x, time, U))
  }
})

test_that("BMP distinguishes the zero-survival convention from Theorem 7", {
  x <- MarkovChain(c(1, 0), matrix(c(0, 1, 1, 0), 2))
  expect_equal(bmp(x, c(0, 1, 2, 5), 1), c(0, 1, 0, 0))
  expect_identical(vcov_bmp(x, 1, 1), 0)
  expect_error(vcov_bmp(x, 2, 1), "positive survival")
  down <- MarkovChain(c(0, 1), x$P)
  expect_equal(bmp(down, 0:2, 1), c(1, 0, 0))
  expect_error(vcov_bmp(down, 1, 1), "positive survival")
  expect_equal(bmp(x, 0:2, integer()), c(1, 0, 0))
  expect_error(vcov_bmp(x, 1, integer()), "positive survival")
  expect_identical(vcov_bmp(x, 0, integer()), 0)
  expect_equal(bmp(x, 0:2, 1:2), rep(0, 3))
  expect_identical(vcov_bmp(x, 2, 1:2), 0)
  single <- MarkovChain(1, matrix(1, 1, 1))
  expect_equal(bmp(single, 0:2, 1), rep(0, 3))
  expect_identical(vcov_bmp(single, 2, 1), 0)
})

test_that("BMP validates inputs and variance irreducibility", {
  x <- MarkovChain(c(0.4, 0.6), matrix(c(0.8, 0.2, 0.3, 0.7), 2,
                                    byrow = TRUE))
  for (f in list(bmp, vcov_bmp)) {
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
  expect_error(vcov_bmp(x, 1:2, 1), "single")
  reducible <- MarkovChain(c(1, 0), diag(2))
  expect_identical(bmp(reducible, 3, 1), 0)
  expect_error(vcov_bmp(reducible, 1, 1), "irreducible")
})
