test_that("sim() returns the requested number of observations", {
  chain <- MarkovChain(a = c("A" = 1, "B" = 0),
                       P = matrix(
                         c(0, 1, 1, 0),
                         nrow = 2,
                         byrow = TRUE,
                         dimnames = list(c("A", "B"), c("A", "B"))
                       ))
  observed <- sim(6, chain)
  expect_length(observed, 6)
})

test_that("sim() follows deterministic transitions exactly", {
  chain <- MarkovChain(a = c("A" = 1, "B" = 0),
                       P = matrix(
                         c(0, 1, 1, 0),
                         nrow = 2,
                         byrow = TRUE,
                         dimnames = list(c("A", "B"), c("A", "B"))
                       ))
  observed <- sim(6, chain)
  expected <- ordered(c("A", "B", "A", "B", "A", "B"), levels = c("A", "B"))
  expect_identical(observed, expected)
})

test_that("sim() draws the first state from the initial distribution", {
  chain <- MarkovChain(a = c("A" = 0, "B" = 1),
                       P = matrix(
                         c(1, 0, 0, 1),
                         nrow = 2,
                         byrow = TRUE,
                         dimnames = list(c("A", "B"), c("A", "B"))
                       ))
  observed <- sim(5, chain)
  expected <- ordered(rep("B", 5), levels = c("A", "B"))
  expect_identical(observed, expected)
})

test_that("sim() returns an ordered factor with all state levels", {
  chain <- MarkovChain(
    a = c("unvisited" = 0, "visited" = 1),
    P = matrix(
      c(1, 0, 0, 1),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("unvisited", "visited"), c("unvisited", "visited"))
    )
  )
  observed <- sim(4, chain)
  expect_s3_class(observed, "factor")
  expect_true(is.ordered(observed))
  expect_identical(levels(observed), c("unvisited", "visited"))
  expect_identical(as.character(observed), rep("visited", 4))
})

test_that("sim() preserves the Markov chain state ordering", {
  chain <- MarkovChain(a = c("z" = 1, "a" = 0),
                       P = matrix(
                         c(0, 1, 1, 0),
                         nrow = 2,
                         byrow = TRUE,
                         dimnames = list(c("z", "a"), c("z", "a"))
                       ))
  observed <- sim(4, chain)
  expected <- ordered(c("z", "a", "z", "a"), levels = c("z", "a"))
  expect_identical(observed, expected)
})

test_that("sim() handles a one-state chain", {
  chain <- MarkovChain(a = c("only" = 1),
                       P = matrix(
                         1,
                         nrow = 1,
                         ncol = 1,
                         dimnames = list("only", "only")
                       ))
  observed <- sim(5, chain)
  expected <- ordered(rep("only", 5), levels = "only")
  expect_identical(observed, expected)
})

test_that("sim() handles a trajectory of length one", {
  chain <- MarkovChain(a = c("A" = 1, "B" = 0),
                       P = matrix(
                         c(0, 1, 1, 0),
                         nrow = 2,
                         byrow = TRUE,
                         dimnames = list(c("A", "B"), c("A", "B"))
                       ))
  observed <- sim(1, chain)
  expected <- ordered("A", levels = c("A", "B"))
  expect_identical(observed, expected)
})

test_that("sim() is reproducible when the random seed is fixed", {
  chain <- MarkovChain(
    a = c("A" = 0.4, "B" = 0.6),
    P = matrix(
      c(0.7, 0.3, 0.2, 0.8),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("A", "B"), c("A", "B"))
    )
  )
  set.seed(123)
  first_path <- sim(20, chain)
  set.seed(123)
  second_path <- sim(20, chain)
  expect_identical(first_path, second_path)
})

test_that("sim() produces only valid states", {
  chain <- MarkovChain(
    a = c("A" = 0.4, "B" = 0.6),
    P = matrix(
      c(0.7, 0.3, 0.2, 0.8),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("A", "B"), c("A", "B"))
    )
  )
  set.seed(321)
  observed <- sim(100, chain)
  expect_false(anyNA(observed))
  expect_true(all(as.character(observed) %in% c("A", "B")))
})

test_that("sim() uses corrected probabilities for tolerated negative roundoff", {
  chain <- MarkovChain(c(1, 0), matrix(c(0, 1, 1, 0), 2, byrow = TRUE))
  chain$P[1, ] <- c(-1e-17, 1)
  observed <- sim(4, chain)
  expect_identical(observed, ordered(c("1", "2", "1", "2"),
                                     levels = c("1", "2")))
  expect_identical(as.numeric(chain$P[1, ]), c(-1e-17, 1))
})

test_that("sim() rejects objects that are not Markov chains", {
  expect_error(sim(5, list(a = 1, P = matrix(1))))
  expect_error(sim(5, NULL))
})

test_that("sim() requires m to be a positive integer", {
  chain <- MarkovChain(a = c("A" = 1, "B" = 0),
                       P = matrix(
                         c(0, 1, 1, 0),
                         nrow = 2,
                         byrow = TRUE,
                         dimnames = list(c("A", "B"), c("A", "B"))
                       ))
  invalid_m <- list(
    0, -1, 1.5, NA_real_, Inf, "5", TRUE, numeric(), c(2, 3), matrix(2)
  )
  for (value in invalid_m) {
    expect_error(sim(value, chain),
                 "`m` must be a single positive integer",
                 info = paste("`m` should be rejected:",
                              paste(value, collapse = ", ")))
  }
})

test_that("sim() returns independent paths in rows in seeded draw order", {
  states <- c("z", "a")
  P <- matrix(c(0.7, 0.3, 0.2, 0.8), 2, byrow = TRUE,
              dimnames = list(states, states))
  chain <- MarkovChain(c(z = 0.4, a = 0.6), P)

  set.seed(123)
  paths <- sim(m = 7, x = chain, n = 4)
  set.seed(123)
  expect_identical(paths, t(replicate(4, as.character(sim(7, chain)))))
  expect_identical(dim(paths), c(4L, 7L))

  set.seed(123)
  path <- sim(7, chain, n = 1)
  expect_identical(path, ordered(paths[1L, ], levels = states))
})

test_that("sim() returns paths that fit() can pool without crossing boundaries", {
  P <- matrix(c(0, 1, 1, 0), 2, byrow = TRUE,
              dimnames = list(c("A", "B"), c("A", "B")))
  chain <- MarkovChain(c(A = 1, B = 0), P)
  paths <- sim(3, chain, n = 5)

  expect_identical(paths, matrix(rep(c("A", "B", "A"), each = 5), 5, 3))
  expect_identical(fit(paths), chain)
  expect_identical(sim(1, chain, n = 5), matrix("A", nrow = 5, ncol = 1))

  single_state <- MarkovChain(c(only = 1), matrix(1, dimnames = list("only", "only")))
  expect_identical(sim(4, single_state, n = 3), matrix("only", 3, 4))
})

test_that("sim() requires n to be a single positive integer", {
  chain <- MarkovChain(c(1, 0), matrix(c(0, 1, 1, 0), 2, byrow = TRUE))
  invalid <- list(0, -1, 1.5, NA_real_, NaN, Inf, "5", TRUE, NULL,
                  numeric(), c(2, 3), matrix(2), .Machine$integer.max + 1)
  for (n in invalid) {
    expect_error(sim(5, chain, n = n), "`n` must be a single positive integer")
  }
})
