test_that("mc_estim() rejects invalid observation sequences", {
  expect_error(
    mc_estim(1),
    "length >= 2"
  )
  expect_error(
    mc_estim(c("a", NA_character_)),
    "contains NA values"
  )
})


test_that("mc_estim() requires exactly one valid init value", {
  expect_error(
    mc_estim(c("a", "b"), init = "invalid")
  )
  expect_error(
    mc_estim(
      c("a", "b"),
      init = c("realized", "uniform")
    )
  )
})


test_that("mc_estim() estimates numeric transitions and sorts states", {
  x <- c(2, 1, 2, 2, 1)
  expected_P <- matrix(c(0, 1, 2/3, 1/3), nrow = 2, byrow = TRUE,
                       dimnames = list(c("1", "2"), c("1", "2")))
  expected_a <- c("1" = 0, "2" = 1)
  expected <- MarkovChain(a = expected_a, P = expected_P)
  observed <- mc_estim(x)
  expect_equal(observed, expected, tolerance = 1e-12)
})


test_that("mc_estim() preserves first-appearance order for character states", {
  x <- c("b", "a", "b", "a")
  expected_P <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE,
                       dimnames = list(c("b", "a"), c("b", "a")))
  expected_a <- c("b" = 0.5, "a" = 0.5)
  expected <- MarkovChain(a = expected_a, P = expected_P)
  observed <- mc_estim(x, init = "uniform")
  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("mc_estim() preserves first-appearance order for unordered factors", {
  x <- factor(c("b", "a", "b", "a"), levels = c("a", "b"))
  expected_P <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE,
                       dimnames = list(c("b", "a"), c("b", "a")))
  expected_a <- c("b" = 1, "a" = 0)
  expected <- MarkovChain(a = expected_a, P = expected_P)
  observed <- mc_estim(x)
  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("mc_estim() sorts ordered-factor states", {
  x <- ordered(c("medium", "low", "high", "medium"),
               levels = c("low", "medium", "high"))
  expected_P <- matrix(
    c(0, 0, 1, 1, 0, 0, 0, 1, 0), nrow = 3, byrow = TRUE,
    dimnames = list(c("low", "medium", "high"), c("low", "medium", "high")))
  expected_a <- c("low" = 0, "medium" = 1, "high"   = 0)
  expected <- MarkovChain(a = expected_a, P = expected_P)
  observed <- mc_estim(x)
  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("mc_estim() treats states with no outgoing observations as absorbing", {
  x <- c("a", "a", "b")
  expected_P <- matrix(c(0.5, 0.5, 0, 1), nrow = 2, byrow = TRUE,
                       dimnames = list(c("a", "b"), c("a", "b")))
  expected_a <- c("a" = 1, "b" = 0)
  expected <- MarkovChain(a = expected_a, P = expected_P)

  expect_warning(observed <- mc_estim(x), "absorbing states: b")
  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("mc_estim() estimates correctly the stationary distribution", {
  # The empirical transition matrix is:
  #       a   b
  #   a   0   1
  #   b 2/3 1/3
  # Its stationary distribution is (0.4, 0.6).
  x <- c("a", "b", "a", "b", "b", "a")
  expected_P <- matrix(c(0, 1, 2/3, 1/3), nrow = 2, byrow = TRUE,
                       dimnames = list(c("a", "b"), c("a", "b")))
  expected_a <- c("a" = 0.4, "b" = 0.6)
  expected <- MarkovChain(a = expected_a, P = expected_P)
  observed <- mc_estim(x, init = "stationary")
  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("all initialization choices agree for a one-state chain", {
  x <- rep("a", 4)
  expected <- MarkovChain(
    a = c("a" = 1),
    P = matrix(1, nrow = 1, ncol = 1, dimnames = list("a", "a")))
  expect_equal(mc_estim(x, init = "realized"), expected)
  expect_equal(mc_estim(x, init = "uniform"), expected)
  expect_equal(mc_estim(x, init = "stationary"), expected)
})
