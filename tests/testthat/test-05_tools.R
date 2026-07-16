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
  expect_error(
    set_state_names(c(0.4, 0.6)),
    "^`x` must be a `StochVec` or `StochMat` object\\.$"
  )

  expect_error(
    set_state_names(matrix(1, nrow = 1)),
    "^`x` must be a `StochVec` or `StochMat` object\\.$"
  )

  expect_error(
    set_state_names(list()),
    "^`x` must be a `StochVec` or `StochMat` object\\.$"
  )
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
