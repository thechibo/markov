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
