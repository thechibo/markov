test_that("vcov_pk() agrees with vcov_pk_sl() for k = 1", {
  a <- c(0.5, 0.5, 0)

  P <- matrix(
    c(
      0.3, 0.7, 0.0,
      0.6, 0.1, 0.3,
      0.5, 0.0, 0.5
    ),
    nrow = 3,
    byrow = TRUE
  )

  x <- MarkovChain(a = a, P = P)

  expect_equal(
    vcov_pk(x, k = 1),
    vcov_pk_sl(x, k = 1),
    tolerance = 1e-10
  )
})


test_that("vcov_pk() agrees with vcov_pk_sl() for several transition horizons", {
  a <- c(0.5, 0.5, 0)

  P <- matrix(
    c(
      0.3, 0.7, 0.0,
      0.6, 0.1, 0.3,
      0.5, 0.0, 0.5
    ),
    nrow = 3,
    byrow = TRUE
  )

  x <- MarkovChain(a = a, P = P)

  for (k in c(2L, 3L, 5L, 10L)) {
    expect_equal(
      vcov_pk(x, k = k),
      vcov_pk_sl(x, k = k),
      tolerance = 1e-10,
      info = paste("Results differ for k =", k)
    )
  }
})


test_that("vcov_pk() and vcov_pk_sl() return matrices of the same dimensions", {
  a <- c(0.5, 0.5, 0)

  P <- matrix(
    c(
      0.3, 0.7, 0.0,
      0.6, 0.1, 0.3,
      0.5, 0.0, 0.5
    ),
    nrow = 3,
    byrow = TRUE
  )

  x <- MarkovChain(a = a, P = P)

  result <- vcov_pk(x, k = 4)
  reference <- vcov_pk_sl(x, k = 4)

  expect_true(is.matrix(result))
  expect_true(is.matrix(reference))

  expect_equal(dim(result), c(9L, 9L))
  expect_identical(dim(result), dim(reference))
  expect_equal(result, reference, tolerance = 1e-10)
})


test_that("vcov_pk() agrees with vcov_pk_sl() for a two-state chain", {
  a <- c(1, 0)

  P <- matrix(
    c(
      0.8, 0.2,
      0.3, 0.7
    ),
    nrow = 2,
    byrow = TRUE
  )

  x <- MarkovChain(a = a, P = P)

  for (k in c(1L, 2L, 4L, 8L)) {
    expect_equal(
      vcov_pk(x, k = k),
      vcov_pk_sl(x, k = k),
      tolerance = 1e-10,
      info = paste("Results differ for k =", k)
    )
  }
})
