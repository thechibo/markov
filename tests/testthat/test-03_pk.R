test_that("k-step covariance functions require a scalar positive integer", {
  chain <- MarkovChain(c(1, 0), matrix(c(0.8, 0.2, 0.3, 0.7),
                                      nrow = 2, byrow = TRUE))
  invalid_k <- list(2.5, 0, -1, NA_real_, NaN, Inf, TRUE, "2",
                    numeric(), c(1, 2), matrix(2), array(2, 1),
                    .Machine$integer.max + 1)
  for (x in list(chain, chain$P)) {
    for (k in invalid_k) {
      expect_error(vcov_pk(x, k), "single positive integer")
      expect_error(vcov_pk_sl(x, k), "single positive integer")
    }
  }
})

test_that("the reference covariance accepts both documented input classes", {
  chain <- MarkovChain(c(1, 0), matrix(c(0.8, 0.2, 0.3, 0.7),
                                      nrow = 2, byrow = TRUE))
  for (k in c(1, 2, 4)) {
    expect_equal(vcov_pk_sl(chain$P, k), vcov_pk_sl(chain, k))
    expect_equal(vcov_pk_sl(chain$P, k), vcov_pk(chain$P, k),
                 tolerance = 1e-12)
    expect_equal(vcov_pk(chain$P, as.integer(k)), vcov_pk(chain$P, k))
  }
  expect_error(vcov_pk_sl(unclass(chain$P), 2), "StochMat")
})

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
