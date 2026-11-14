test_that("maintainability includes recovery by time and initial up mass", {
  P <- matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE)
  x <- MarkovChain(c(0.4, 0.6), P)
  time <- c(3, 0, 1, 2, 1)
  expect_equal(main(x, time, down = 2), 1 - 0.6 * 0.7^time)
  expect_gt(main(x, 2, 2), avail(x, 2, 1))
  starting_down <- MarkovChain(c(0, 1), P)
  expect_equal(main(starting_down, 0:3, 2), 1 - 0.7^(0:3))
  for (f in list(vcov_main, vcov_main_sl)) {
    for (t in c(1, 2, 5)) {
      expected <- (0.6 * t * 0.7^(t - 1))^2 * 0.7 * 0.3 / 0.4
      expect_equal(f(x, t, 2), expected)
    }
    expect_identical(f(x, 0, 2), 0)
  }
})

test_that("maintainability preserves noncontiguous down-state selection", {
  P <- rbind(c(0.6, 0.3, 0.1), c(0.2, 0.5, 0.3), c(0.4, 0.1, 0.5))
  labels <- c("failed", "working", "repair")
  dimnames(P) <- list(labels, labels)
  a <- setNames(c(0.7, 0.2, 0.1), labels)
  x <- MarkovChain(a, P)
  down <- c("repair", "failed", "repair")
  Q <- P[c(1, 3), c(1, 3)]
  expect_equal(main(x, 3, down), 1 - sum(a[c(1, 3)] %*% Q %*% Q %*% Q))
  expect_equal(main(x, 0:3, down) + dssc(x, 0:3, c(1, 3)), rep(1, 4))
  expect_equal(vcov_main(x, 3, down), vcov_dssc(x, 3, c(1, 3)))
  expect_equal(vcov_main_sl(x, 3, down), vcov_dssc_sl(x, 3, c(1, 3)))
})

test_that("maintainability preserves boundary cases and validation", {
  x <- MarkovChain(c(0.4, 0.6), matrix(c(0.8, 0.2, 0.3, 0.7), 2,
                                    byrow = TRUE))
  expect_equal(main(x, 0:2, integer()), rep(1, 3))
  expect_equal(main(x, 0:2, 1:2), rep(0, 3))
  reducible <- MarkovChain(c(0, 1), diag(2))
  expect_equal(main(reducible, 0:2, 2), rep(0, 3))
  for (f in list(vcov_main, vcov_main_sl)) {
    expect_identical(f(x, 3, integer()), 0)
    expect_identical(f(x, 3, 1:2), 0)
    expect_error(f(x, 1:2, 2), "single")
    expect_error(f(reducible, 1, 2), "irreducible")
  }
  for (f in list(main, vcov_main, vcov_main_sl)) {
    expect_error(f(x, NA_real_, 2), "time")
    expect_error(f(x, -1, 2), "time")
    expect_error(f(x, 1, "unknown"), "unknown state")
    expect_error(f(x, 1, 2.5), "indices")
    invalid <- x
    invalid$a[] <- c(2, 0)
    expect_error(f(invalid, 1, 2), "valid StochVec")
  }
})
