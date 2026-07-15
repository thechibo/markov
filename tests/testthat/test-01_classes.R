# Tests for probability and Markov-chain classes ---------------------------

make_test_matrix <- function() {
  matrix(c(0.7, 0.3,
           0.2, 0.8), nrow = 2, byrow = TRUE,
         dimnames = list(c("A", "B"), c("A", "B")))
}

make_test_chain <- function(a = c("A" = 0.6, "B" = 0.4)) {
  MarkovChain(a = a, P = make_test_matrix())
}


# StochVec -----------------------------------------------------------------

test_that("StochVec() constructs a named stochastic vector", {
  observed <- StochVec(c("A" = 0.25, "B" = 0.75))

  expect_s3_class(observed, "StochVec")
  expect_true(is.numeric(observed))
  expect_equal(as.numeric(observed), c(0.25, 0.75))
  expect_identical(names(observed), c("A", "B"))
  expect_true(is.StochVec(observed))
})

test_that("StochVec() accepts zero probabilities", {
  observed <- StochVec(c("A" = 1, "B" = 0))

  expect_s3_class(observed, "StochVec")
  expect_equal(as.numeric(observed), c(1, 0))
})

test_that("StochVec() rejects missing and negative values", {
  expect_error(StochVec(c(NA_real_, 1)), "non-missing")
  expect_error(StochVec(c(-0.1, 1.1)), "greater than zero")
})

test_that("StochVec() rejects vectors that do not sum to one", {
  expect_error(StochVec(c(0.2, 0.7)), "sum to 1")
  expect_error(StochVec(c(0.5, 0.6)), "sum to 1")
})

test_that("StochVec() preserves an existing valid object", {
  x <- StochVec(c("A" = 0.4, "B" = 0.6))

  expect_equal(StochVec(x), x)
  expect_true(is.StochVec(x))
  expect_false(is.StochVec(as.numeric(x)))
})

test_that("print.StochVec() prints a summary and returns invisibly", {
  x <- StochVec(c("A" = 0.4, "B" = 0.6))
  output <- capture.output(returned <- print(x))

  expect_match(paste(output, collapse = "\n"), "Stochastic Vector")
  expect_match(paste(output, collapse = "\n"), "length: 2")
  expect_identical(returned, x)
})


# StochMat -----------------------------------------------------------------

test_that("StochMat() constructs a stochastic matrix", {
  observed <- StochMat(make_test_matrix())

  expect_s3_class(observed, "StochMat")
  expect_true(is.matrix(observed))
  expect_true(is.numeric(observed))
  expect_identical(dim(observed), c(2L, 2L))
  expect_equal(rowSums(observed), c(1, 1), tolerance = 1e-12)
  expect_identical(rownames(observed), c("A", "B"))
  expect_identical(colnames(observed), c("A", "B"))
  expect_true(is.StochMat(observed))
})

test_that("StochMat() can construct a matrix from a vector", {
  observed <- StochMat(c(0.7, 0.3, 0.2, 0.8), nrow = 2, byrow = TRUE)

  expect_s3_class(observed, "StochMat")
  expect_identical(dim(observed), c(2L, 2L))
  expect_equal(rowSums(observed), c(1, 1), tolerance = 1e-12)
})

test_that("StochMat() accepts zero transition probabilities", {
  observed <- StochMat(diag(2))

  expect_s3_class(observed, "StochMat")
  expect_equal(unclass(observed), diag(2))
})

test_that("StochMat() rejects nonsquare and nonnumeric matrices", {
  expect_error(StochMat(matrix(1 / 3, nrow = 2, ncol = 3)), "square matrix")
  expect_error(StochMat(matrix(letters[1:4], nrow = 2)), "numeric matrix")
})

test_that("StochMat() rejects missing and nonfinite values", {
  P_na <- make_test_matrix()
  P_na[1, 1] <- NA_real_

  P_inf <- make_test_matrix()
  P_inf[1, 1] <- Inf

  expect_error(StochMat(P_na), "finite and non-missing")
  expect_error(StochMat(P_inf), "finite and non-missing")
})

test_that("StochMat() rejects materially negative values", {
  P <- matrix(c(-0.01, 1.01,
                0.20, 0.80), nrow = 2, byrow = TRUE)

  expect_error(StochMat(P), "must be >= 0")
})

test_that("StochMat() clips tiny negative values and renormalizes rows", {
  P <- matrix(c(-1e-13, 1 + 1e-13,
                0.4, 0.6), nrow = 2, byrow = TRUE)

  observed <- StochMat(P)

  expect_equal(observed[1, ], c(0, 1), tolerance = 1e-12)
  expect_equal(rowSums(observed), c(1, 1), tolerance = 1e-12)
})

test_that("StochMat() rejects rows that do not sum to one", {
  P <- matrix(c(0.7, 0.2,
                0.2, 0.8), nrow = 2, byrow = TRUE)

  expect_error(StochMat(P), "Row values must sum to 1")
})

test_that("is.StochMat() distinguishes matrices from StochMat objects", {
  P <- make_test_matrix()

  expect_false(is.StochMat(P))
  expect_true(is.StochMat(StochMat(P)))
})

test_that("print.StochMat() prints a summary and returns invisibly", {
  P <- StochMat(make_test_matrix())
  output <- capture.output(returned <- print(P))

  expect_match(paste(output, collapse = "\n"), "Stochastic Matrix")
  expect_match(paste(output, collapse = "\n"), "dim: 2 x 2")
  expect_identical(returned, P)
})


# MarkovChain --------------------------------------------------------------

test_that("MarkovChain() constructs its component classes", {
  observed <- make_test_chain()

  expect_s3_class(observed, "MarkovChain")
  expect_s3_class(observed$a, "StochVec")
  expect_s3_class(observed$P, "StochMat")
  expect_identical(names(observed), c("a", "P"))
  expect_true(is.MarkovChain(observed))
})

test_that("MarkovChain() accepts zero initial and transition probabilities", {
  observed <- MarkovChain(a = c("A" = 1, "B" = 0), P = diag(2))

  expect_equal(as.numeric(observed$a), c(1, 0))
  expect_equal(unclass(observed$P), diag(2))
})

test_that("MarkovChain() propagates component validation errors", {
  expect_error(MarkovChain(a = c(0.4, 0.4), P = make_test_matrix()), "sum to 1")

  invalid_P <- matrix(c(0.7, 0.1,
                        0.2, 0.8), nrow = 2, byrow = TRUE)
  expect_error(MarkovChain(a = c(0.5, 0.5), P = invalid_P),
               "Row values must sum to 1")
})

test_that("is.MarkovChain() distinguishes lists from MarkovChain objects", {
  chain <- make_test_chain()
  plain <- unclass(chain)

  expect_true(is.MarkovChain(chain))
  expect_false(is.MarkovChain(plain))
})

test_that("print.MarkovChain() prints both components and returns invisibly", {
  chain <- make_test_chain()
  output <- capture.output(returned <- print(chain))
  text <- paste(output, collapse = "\n")

  expect_match(text, "Markov Chain")
  expect_match(text, "Number of states: 2")
  expect_match(text, "Initial distribution")
  expect_match(text, "Transition matrix")
  expect_identical(returned, chain)
})


# MarkovChainList ----------------------------------------------------------

test_that("MarkovChainList() constructs a list of Markov chains", {
  first <- make_test_chain(c("A" = 1, "B" = 0))
  second <- make_test_chain(c("A" = 0, "B" = 1))
  observed <- MarkovChainList(list(first = first, second = second))

  expect_s3_class(observed, "MarkovChainList")
  expect_length(observed, 2)
  expect_identical(names(observed), c("first", "second"))
  expect_true(all(vapply(observed, is.MarkovChain, logical(1))))
})

test_that("MarkovChainList() rejects nonlists and invalid elements", {
  expect_error(MarkovChainList(make_test_chain()))
  expect_error(MarkovChainList(list(make_test_chain(), "not a chain")),
               "MarkovChain")
})

test_that("subsetting MarkovChainList preserves its class", {
  chains <- MarkovChainList(list(
    first = make_test_chain(c("A" = 1, "B" = 0)),
    second = make_test_chain(c("A" = 0, "B" = 1))
  ))

  observed <- chains[2]

  expect_s3_class(observed, "MarkovChainList")
  expect_length(observed, 1)
  expect_identical(names(observed), "second")
  expect_s3_class(observed[[1]], "MarkovChain")
})

test_that("MarkovChainList() accepts an empty list", {
  observed <- MarkovChainList(list())

  expect_s3_class(observed, "MarkovChainList")
  expect_length(observed, 0)
})


# Invariant regression tests -----------------------------------------------
# These tests expose validation gaps in the current implementation.

test_that("MarkovChain requires exactly the elements a and P", {
  a <- StochVec(c("A" = 0.5, "B" = 0.5))
  P <- StochMat(make_test_matrix())

  malformed <- structure(
    list(a = a, P = P, a = a),
    class = c("MarkovChain", "list")
  )

  expect_error(validate_MarkovChain(malformed))
})

test_that("MarkovChain requires compatible vector and matrix dimensions", {
  expect_error(
    MarkovChain(a = c(0.5, 0.5), P = diag(3))
  )
})

test_that("MarkovChainList validates objects rather than class labels only", {
  malformed <- structure(
    list(
      a = StochVec(c(0.5, 0.5)),
      P = StochMat(diag(3))
    ),
    class = c("MarkovChain", "list")
  )

  expect_error(MarkovChainList(list(malformed)))
})
