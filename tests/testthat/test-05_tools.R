# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set_state_names()                ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test_that("set_state_names() assigns sequential names to an unnamed StochVec", {
  x <- StochVec(c(0.4, 0.6))
  result <- set_state_names(x)
  expect_identical(names(result), c("1", "2"))
  expect_equal(as.numeric(result), c(0.4, 0.6))
  expect_identical(class(result), class(x))
})

test_that("set_state_names() preserves existing StochVec names", {
  x <- StochVec(c(low = 0.4, high = 0.6))
  result <- set_state_names(x)
  expect_identical(names(result), c("low", "high"))
  expect_equal(as.numeric(result), c(0.4, 0.6))
  expect_identical(class(result), class(x))
})

test_that("set_state_names.StochVec ignores additional arguments", {
  x <- StochVec(c(0.25, 0.75))
  result <- set_state_names(x, unused_argument = TRUE)
  expect_identical(names(result), c("1", "2"))
})

test_that("set_state_names() assigns sequential names to an unnamed StochMat", {
  x <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE))
  result <- set_state_names(x)
  expect_identical(rownames(result), c("1", "2"))
  expect_identical(colnames(result), c("1", "2"))
  expect_equal(unname(result), unname(x))
  expect_identical(class(result), class(x))
})

test_that("set_state_names() preserves identical StochMat dimnames", {
  state_names <- c("low", "high")
  x <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE,
                       dimnames = list(state_names, state_names)))
  result <- set_state_names(x)
  expect_identical(rownames(result), state_names)
  expect_identical(colnames(result), state_names)
  expect_equal(unname(result), unname(x))
})

test_that("set_state_names.StochMat ignores additional arguments", {
  x <- StochMat(diag(2))
  result <- set_state_names(x, unused_argument = TRUE)
  expect_identical(rownames(result), c("1", "2"))
  expect_identical(colnames(result), c("1", "2"))
})

test_that("set_state_names() rejects unsupported objects", {
  expect_error(set_state_names(c(0.4, 0.6)),
    "^`x` must be a `StochVec` or `StochMat` object\\.$")
  expect_error(set_state_names(matrix(1, nrow = 1)),
    "^`x` must be a `StochVec` or `StochMat` object\\.$")
  expect_error(set_state_names(list()),
    "^`x` must be a `StochVec` or `StochMat` object\\.$")
})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get_state_names()                ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test_that("get_state_names() returns StochVec names", {
  x <- StochVec(c(low = 0.4, high = 0.6))
  expect_identical(get_state_names(x), c("low", "high"))
})

test_that("get_state_names() returns sequential names for an unnamed StochVec", {
  x <- StochVec(c(0.4, 0.6))
  expect_identical(get_state_names(x), c("1", "2"))
})

test_that("get_state_names.StochVec ignores additional arguments", {
  x <- StochVec(c(low = 0.4, high = 0.6))
  expect_identical(get_state_names(x, unused_argument = TRUE), c("low", "high"))
})

test_that("get_state_names() returns StochMat row names", {
  x <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE,
                       dimnames = list(c("low", "high"), c("low", "high"))))
  expect_identical(get_state_names(x), c("low", "high"))
})

test_that("get_state_names() returns sequential names for an unnamed StochMat", {
  x <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE))
  expect_identical(get_state_names(x), c("1", "2"))
})

test_that("get_state_names() returns state names for a MarkovChain", {
  a <- StochVec(c("low" = 0.3, "high" = 0.7))
  P <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE,
                       dimnames = list(c("low", "high"), c("low", "high"))))
  x <- MarkovChain(a, P)
  expect_identical(get_state_names(x), c("low", "high"))
})

test_that("get_state_names.MarkovChain ignores additional arguments", {
  a <- StochVec(c("low" = 0.3, "high" = 0.7))
  P <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE,
                       dimnames = list(c("low", "high"), c("low", "high"))))
  x <- MarkovChain(a, P)
  expect_identical(get_state_names(x, unused_argument = TRUE), c("low", "high"))
})

test_that("get_state_names() rejects unsupported objects", {
  expected_message <- paste0(
    "^`x` must be a `StochVec`, a `StochMat`, ",
    "or a `MarkovChain` object\\.$"
  )
  expect_error(get_state_names(c(low = 0.4, high = 0.6)), expected_message)
  expect_error(get_state_names(matrix(1, nrow = 1)), expected_message)
  expect_error(get_state_names(list()), expected_message)
})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# nstates()                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test_that("nstates() returns the length of a StochVec", {
  a <- StochVec(c(low = 0.4, high = 0.6))
  expect_identical(nstates(a), 2L)
})

test_that("nstates() returns the number of rows in a StochMat", {
  P <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE))
  expect_identical(nstates(P), 2L)
})

test_that("nstates() returns the number of rows in a MarkovChain", {
  a <- StochVec(c(0.4, 0.6))
  P <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE))
  x <- MarkovChain(a, P)
  expect_identical(nstates(x), 2L)
})

test_that("nstates() rejects unsupported objects", {
  expect_error(nstates(c(0.4, 0.6)),
    "`x` must be a `StochVec`, a `StochMat`, or a `MarkovChain` object.")
  expect_error(nstates(matrix(1, nrow = 1)),
    "`x` must be a `StochVec`, a `StochMat`, or a `MarkovChain` object.")
})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Classes                          ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

make_closed_classes_matrix <- function() {
  states <- c("a", "b", "c")
  StochMat(matrix(
    c(0.8, 0.2, 0,
      0.3, 0.7, 0,
      0.0, 0.0, 1),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(states, states)))
}

make_transient_class_matrix <- function() {
  states <- c("a", "b", "c")
  StochMat(matrix(
      c(0.5, 0.4, 0.1,
        0.4, 0.5, 0.1,
        0.0, 0.0, 1.0),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(states, states)))
}

test_that("mc_classes() returns communicating classes using state names", {
  P <- make_closed_classes_matrix()
  expect_identical(mc_classes(P), list(c("a", "b"), "c"))
})

test_that("mc_classes() dispatches on MarkovChain objects", {
  P <- make_closed_classes_matrix()
  a <- StochVec(c(a = 1, b = 0, c = 0))
  x <- MarkovChain(a = a, P = P)
  expect_identical(mc_classes(x), mc_classes(P))
})

test_that("graph methods validate transition matrices and Markov chains", {
  P <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE))
  chain <- MarkovChain(c(0.6, 0.4), P)
  chain$a <- NULL
  P[1L, 1L] <- 0.5

  for (graph_method in list(mc_classes, is_irreducible)) {
    expect_error(graph_method(chain), "exactly two elements")
    expect_error(graph_method(P), "Row values must sum to 1")
  }
})

test_that("Ppowers() validates inputs and preserves state names", {
  P <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE,
                       dimnames = list(c("a", "b"), c("a", "b"))))
  powers <- Ppowers(P, 2L)

  expect_length(powers, 3L)
  expect_identical(dimnames(powers[[1L]]), dimnames(P))
  expect_equal(powers[[3L]], P %*% P)
  expect_error(Ppowers(matrix(1:6, 2L), 2L), "square numeric matrix")
  expect_error(Ppowers(P, 2.5), "single positive integer")
})

test_that("is_com_class() identifies communicating classes", {
  P <- make_closed_classes_matrix()
  expect_true(is_com_class(P, c("a", "b")))
  expect_true(is_com_class(P, "c"))
  expect_false(is_com_class(P, c("a", "c")))
  expect_false(is_com_class(P, "a"))
})

test_that("is_com_class() accepts integer values for default state names", {
  P <- StochMat(matrix(
      c(0.5, 0.5, 0,
        0.5, 0.5, 0,
        0.0, 0.0, 1),
      nrow = 3,
      byrow = TRUE))
  expect_true(is_com_class(P, c(1L, 2L)))
  expect_true(is_com_class(P, 3L))
  expect_false(is_com_class(P, c(1L, 3L)))
})

test_that("is_com_class() validates its inputs", {
  P <- make_closed_classes_matrix()
  expect_error(is_com_class(diag(3), "a"),
               "`P` must be a `StochMat` object", fixed = TRUE)
  expect_error(is_com_class(P, character()),
               "`c` must be a non-empty vector", fixed = TRUE)
  expect_error(is_com_class(P, NA_character_),
    "`c` must not contain missing values|`c` contains NAs")
  expect_error(is_com_class(P, "unknown"),
    "All values in `c` must be valid state names", fixed = TRUE)
})

test_that("is_rec_class() distinguishes recurrent and transient classes", {
  P <- make_transient_class_matrix()
  expect_false(is_rec_class(P, c("a", "b")))
  expect_true(is_rec_class(P, "c"))
})

test_that("is_rec_class() requires a complete communicating class", {
  P <- make_transient_class_matrix()
  expect_error(is_rec_class(P, "a"),
               "`c` must identify a communicating class", fixed = TRUE)
  expect_error(is_rec_class(P, c("a", "c")),
               "`c` must identify a communicating class", fixed = TRUE)
})

test_that("class helpers use the state names assigned during validation", {
  P <- new_StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE))
  states <- c("1", "2")

  expect_true(is_com_class(P, states))
  expect_true(is_rec_class(P, states))
  expect_equal(stat_on_class(P, states), stat(StochMat(P)))
})

test_that("stat_on_class() computes and embeds a stationary distribution", {
  P <- make_closed_classes_matrix()
  result <- stat_on_class(P, c("a", "b"))
  expect_s3_class(result, "StochVec")
  expect_identical(names(result), c("a", "b", "c"))
  expect_equal(as.numeric(result), c(0.6, 0.4, 0), tolerance = 1e-8)
  expect_equal(sum(result), 1, tolerance = 1e-8)
})

test_that("stat_on_class() handles an absorbing state", {
  P <- make_closed_classes_matrix()
  result <- stat_on_class(P, "c")
  expect_s3_class(result, "StochVec")
  expect_identical(names(result), c("a", "b", "c"))
  expect_equal(as.numeric(result), c(0, 0, 1))
})

test_that("stat_on_class() rejects a transient communicating class", {
  P <- make_transient_class_matrix()
  expect_error(stat_on_class(P, c("a", "b")),
               "`c` must identify a recurrent class", fixed = TRUE)
})

test_that("stationary distribution is invariant on its class", {
  P <- make_closed_classes_matrix()
  result <- stat_on_class(P, c("a", "b"))
  expect_equal(as.numeric(result %*% P), as.numeric(result), tolerance = 1e-8)
})
