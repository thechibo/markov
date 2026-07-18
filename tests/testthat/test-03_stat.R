make_stoch_mat <- function(values, states) {
  StochMat(matrix(
    values,
    nrow = length(states),
    byrow = TRUE,
    dimnames = list(states, states)
  ))
}

expect_stationary <- function(pi, P) {
  expect_s3_class(pi, "StochVec")
  expect_equal(sum(pi), 1, tolerance = 1e-10)
  expect_true(all(pi >= -1e-10))
  expect_equal(as.numeric(pi %*% P), as.numeric(pi), tolerance = 1e-10)
}

test_that("stat() computes the unique stationary distribution", {
  P <- make_stoch_mat(c(0.8, 0.2, 0.3, 0.7), c("A", "B"))
  observed <- stat(P)
  expect_stationary(observed, P)
  expect_equal(as.numeric(observed), c(0.6, 0.4), tolerance = 1e-12)
  expect_identical(names(observed), c("A", "B"))
})

test_that("stat() returns a unique distribution for a periodic chain", {
  P <- make_stoch_mat(c(0, 1, 1, 0), c("A", "B"))
  observed <- stat(P)
  expect_stationary(observed, P)
  expect_equal(as.numeric(observed), c(0.5, 0.5), tolerance = 1e-12)
})

test_that("stat() assigns zero probability to transient states", {
  P <- make_stoch_mat(c(0, 1, 0, 0, 0, 1, 0, 1, 0), c("transient", "A", "B"))
  observed <- stat(P)
  expect_stationary(observed, P)
  expect_equal(as.numeric(observed), c(0, 0.5, 0.5), tolerance = 1e-12)
  expect_identical(names(observed), c("transient", "A", "B"))
})

test_that("stat() dispatches correctly for MarkovChain objects", {
  P <- make_stoch_mat(c(0.8, 0.2, 0.3, 0.7), c("A", "B"))
  chain <- MarkovChain(a = c("A" = 1, "B" = 0), P = P)
  expect_equal(stat(chain), stat(P), tolerance = 1e-12)
})

test_that("the stationary distribution does not depend on the initial law", {
  P <- make_stoch_mat(c(0.8, 0.2, 0.3, 0.7), c("A", "B"))
  first <- MarkovChain(a = c("A" = 1, "B" = 0), P = P)
  second <- MarkovChain(a = c("A" = 0, "B" = 1), P = P)
  expect_equal(stat(first), stat(second), tolerance = 1e-12)
})

test_that("stat() returns one extremal distribution per closed class", {
  P <- make_stoch_mat(c(1, 0, 0, 0, 1, 0, 0.4, 0.6, 0), c("A", "B", "transient"))
  observed <- stat(P)
  expect_type(observed, "list")
  expect_length(observed, 2)
  expect_true(all(vapply(observed, inherits, logical(1), what = "StochVec")))

  for (pi in observed) {
    expect_stationary(pi, P)
    expect_equal(as.numeric(pi["transient"]), 0, tolerance = 1e-12)
  }

  supports <- vapply(observed,
                     function(pi) { names(pi)[which.max(pi)] },
                     character(1))
  expect_setequal(supports, c("A", "B"))

  pi_A <- observed[[which(supports == "A")]]
  pi_B <- observed[[which(supports == "B")]]

  expect_equal(as.numeric(pi_A), c(1, 0, 0), tolerance = 1e-12)
  expect_equal(as.numeric(pi_B), c(0, 1, 0), tolerance = 1e-12)
})

test_that("vcov_stat() returns a valid covariance matrix", {
  P <- make_stoch_mat(c(0.8, 0.2, 0.3, 0.7), c("A", "B"))
  chain <- MarkovChain(a = c("A" = 1, "B" = 0), P = P)
  observed <- vcov_stat(chain)
  expect_type(observed, "double")
  expect_identical(dim(observed), c(2L, 2L))
  expect_true(isSymmetric(observed, tol = 1e-10))
  expect_equal(rowSums(observed), c(0, 0), tolerance = 1e-10)
  expect_equal(colSums(observed), c(0, 0), tolerance = 1e-10)
  expect_true(all(
    eigen(observed, symmetric = TRUE, only.values = TRUE)$values >= -1e-10
  ))
})

test_that("vcov_stat() rejects chains with non-unique stationary distributions", {
  P <- make_stoch_mat(c(1, 0, 0, 1), c("A", "B"))
  chain <- MarkovChain(a = c("A" = 0.5, "B" = 0.5), P = P)
  expect_error(vcov_stat(chain))
})

test_that("vcov_stat() accepts a StochMat object as documented", {
  P <- make_stoch_mat(c(0.8, 0.2, 0.3, 0.7), c("A", "B"))
  chain <- MarkovChain(a = c("A" = 1, "B" = 0), P = P)
  expect_equal(vcov_stat(P), vcov_stat(chain), tolerance = 1e-12)
})
