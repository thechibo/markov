test_that("fit() rejects invalid observation sequences", {
  expect_error(
    fit(1),
    "length >= 2"
  )
  expect_error(
    fit(c("a", NA_character_)),
    "contains NA values"
  )
  for (x in list(array(1:4, 4), list("a", "b"))) {
    expect_error(fit(x), "must be an atomic vector")
  }
  for (x in list(c("a", ""), factor(c("", "a")))) {
    expect_error(fit(x), "State labels in `x` must be .*non-empty")
  }
})

test_that("fit pools paths without between-path transitions", {
  x <- rbind(c("b", "a", "a"), c("a", "b", "b"), c("a", "a", "b"))
  fitted <- fit(x)
  expect_equal(as.numeric(fitted$a), c(1 / 3, 2 / 3))
  expect_equal(rownames(fitted$P), c("b", "a"))
  expect_equal(unname(unclass(fitted$P)), matrix(0.5, 2, 2))
  expect_equal(as.numeric(fit(x, "uniform")$a), c(0.5, 0.5))
  expect_equal(as.numeric(fit(x, "stationary")$a), c(0.5, 0.5))
  expect_equal(fit(matrix(x[1, ], nrow = 1)), fit(x[1, ]))
  expect_equal(fit(matrix(1, 3, 4))$a, StochVec(c("1" = 1)))
  expect_warning(terminal <- fit(rbind(c(10, 20), c(10, 20))), "absorbing")
  expect_equal(as.numeric(terminal$a), c(1, 0))
  expect_equal(unname(unclass(terminal$P)), matrix(c(0, 1, 0, 1), 2, byrow = TRUE))
  close <- c(1, 1 + 1e-15)
  expect_equal(nrow(fit(rbind(close, rev(close)))$P), 2L)
})

test_that("trajectory matrices are validated without coercion", {
  x <- matrix(c("a", "b", "b", "a"), 2)
  expect_identical(validate_trajectories(x), x)
  f <- factor(c("a", "b"))
  expect_identical(validate_trajectories(f), f)
  for (bad in list(matrix(numeric(), 0, 3), matrix(1, 3, 1), matrix(1, 2, 0))) {
    expect_error(fit(bad), "at least one row and two columns")
  }
  expect_error(fit(matrix(c(1, NA, 2, 1), 2)), "NA values")
  expect_error(fit(matrix(list(1, 2, 1, 2), 2)), "atomic vector")
  expect_error(fit(array(1, c(2, 2, 2))), "atomic vector")
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
