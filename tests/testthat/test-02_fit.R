test_that("fit() rejects invalid observation sequences", {
  expect_error(
    fit(1),
    "length >= 2"
  )
  expect_error(
    fit(c("a", NA_character_)),
    "contains NA values"
  )
  for (x in list(matrix(1:4, nrow = 2), array(1:4, 4), list("a", "b"))) {
    expect_error(fit(x), "must be an atomic vector")
  }
  for (x in list(c("a", ""), factor(c("", "a")))) {
    expect_error(fit(x), "State labels in `x` must be .*non-empty")
  }
})

test_that("fit() requires exactly one valid init value", {
  invalid <- list(
    "invalid", "r", NULL, NA_character_, c("realized", "uniform"),
    c("realized", "uniform", "stationary")
  )
  for (init in invalid) {
    expect_error(fit(c("a", "b"), init = init), "`init` must be one of")
  }
})

test_that("fit() estimates numeric transitions and sorts states", {
  x <- c(2, 1, 2, 2, 1)
  expected_P <- matrix(c(0, 1, 2/3, 1/3), nrow = 2, byrow = TRUE,
                       dimnames = list(c("1", "2"), c("1", "2")))
  expected_a <- c("1" = 0, "2" = 1)
  expected <- MarkovChain(a = expected_a, P = expected_P)
  observed <- fit(x)
  expect_equal(observed, expected, tolerance = 1e-12)
})


test_that("fit() distinguishes numeric states with colliding printed labels", {
  for (states in list(c(1, 1 + 1e-15), c(1e10, 1e10 + 1e-5))) {
    x <- states[c(2, 1, 2, 1)]
    for (init in c("realized", "uniform", "stationary")) {
      observed <- fit(x, init = init)
      expect_equal(as.numeric(rownames(observed$P)), states, tolerance = 0)
      expect_identical(anyDuplicated(rownames(observed$P)), 0L)
      expect_identical(names(observed$a), rownames(observed$P))
      expect_equal(unname(unclass(observed$P)), matrix(c(0, 1, 1, 0), 2))
      expected_a <- if (init == "realized") c(0, 1) else c(0.5, 0.5)
      expect_equal(as.numeric(observed$a), expected_a)
    }
    expect_warning(terminal <- fit(states[c(1, 1, 2)]), "absorbing states")
    expect_equal(unname(unclass(terminal$P)),
                 matrix(c(0.5, 0.5, 0, 1), 2, byrow = TRUE))
  }
})

test_that("fit() preserves first-appearance order for character states", {
  x <- c("b", "a", "b", "a")
  expected_P <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE,
                       dimnames = list(c("b", "a"), c("b", "a")))
  expected_a <- c("b" = 0.5, "a" = 0.5)
  expected <- MarkovChain(a = expected_a, P = expected_P)
  observed <- fit(x, init = "uniform")
  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("fit() preserves first-appearance order for unordered factors", {
  x <- factor(c("b", "a", "b", "a"), levels = c("a", "b"))
  expected_P <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE,
                       dimnames = list(c("b", "a"), c("b", "a")))
  expected_a <- c("b" = 1, "a" = 0)
  expected <- MarkovChain(a = expected_a, P = expected_P)
  observed <- fit(x)
  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("fit() sorts ordered-factor states", {
  x <- ordered(c("medium", "low", "high", "medium"),
               levels = c("low", "medium", "high"))
  expected_P <- matrix(
    c(0, 0, 1, 1, 0, 0, 0, 1, 0), nrow = 3, byrow = TRUE,
    dimnames = list(c("low", "medium", "high"), c("low", "medium", "high")))
  expected_a <- c("low" = 0, "medium" = 1, "high"   = 0)
  expected <- MarkovChain(a = expected_a, P = expected_P)
  observed <- fit(x)
  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("fit() treats states with no outgoing observations as absorbing", {
  x <- c("a", "a", "b")
  expected_P <- matrix(c(0.5, 0.5, 0, 1), nrow = 2, byrow = TRUE,
                       dimnames = list(c("a", "b"), c("a", "b")))
  expected_a <- c("a" = 1, "b" = 0)
  expected <- MarkovChain(a = expected_a, P = expected_P)

  expect_warning(observed <- fit(x), "absorbing states: b")
  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("fit() makes terminal-only numeric states absorbing by position", {
  for (terminal in c(20, 0, -5, 0.5)) {
    labels <- as.character(c(10, terminal))
    expected_P <- matrix(c(0.5, 0.5, 0, 1), nrow = 2, byrow = TRUE,
                         dimnames = list(labels, labels))
    expected_a <- setNames(c(1, 0), labels)
    state_order <- order(c(10, terminal))
    expected <- MarkovChain(
      a = expected_a[state_order],
      P = expected_P[state_order, state_order]
    )

    expect_warning(
      observed <- fit(c(10, 10, terminal)),
      paste("absorbing states:", terminal),
      fixed = TRUE
    )
    expect_equal(observed, expected, tolerance = 1e-12,
                 info = paste("Terminal state:", terminal))
  }
})

test_that("fit() estimates correctly the stationary distribution", {
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
  observed <- fit(x, init = "stationary")
  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("all initialization choices agree for a one-state chain", {
  x <- rep("a", 4)
  expected <- MarkovChain(
    a = c("a" = 1),
    P = matrix(1, nrow = 1, ncol = 1, dimnames = list("a", "a")))
  expect_equal(fit(x, init = "realized"), expected)
  expect_equal(fit(x, init = "uniform"), expected)
  expect_equal(fit(x, init = "stationary"), expected)
})
