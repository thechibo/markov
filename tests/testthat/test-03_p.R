test_that("covariance blocks retain rare transition probabilities", {
  rare <- 1e-18
  expected <- matrix(c(rare, -rare, -rare, rare), nrow = 2)

  expect_identical(lambda_block(c(1, rare)), expected)
  expect_identical(lambda_block(c(rare, 1)), expected)
  expect_equal(lambda_block(1), matrix(0, nrow = 1))
})

test_that("covariance functions reject corrupted classed inputs", {
  chain <- MarkovChain(rep(1 / 3, 3), matrix(
    c(0.8, 0.1, 0.1, 0.2, 0.5, 0.3, 0.3, 0.3, 0.4), 3, byrow = TRUE
  ))
  covariance_functions <- list(
    vcov_p, function(x) vcov_pk(x, 2), vcov_stat,
    function(x) vcov_pk_sl(x, 2)
  )
  for (row in list(c(0.8, 0.3, -0.1), c(0.8, 0.3, 0.1),
                   c(NA_real_, 0.1, 0.1), c(Inf, 0.1, 0.1))) {
    invalid <- chain
    invalid$P[1, ] <- row
    for (covariance in covariance_functions) {
      expect_error(covariance(invalid$P))
      expect_error(covariance(invalid), "valid StochMat")
    }
  }
  invalid <- chain
  invalid$a[1] <- -1
  for (covariance in covariance_functions) {
    expect_error(covariance(invalid), "valid StochVec")
  }
})

test_that("covariance functions use the cleaned transition matrix", {
  P <- StochMat(matrix(c(0, 0.4, 0.6, 0.2, 0.5, 0.3, 0.3, 0.3, 0.4),
                       3, byrow = TRUE))
  roundoff <- P
  roundoff[1, 1] <- -1e-17
  roundoff[2, ] <- roundoff[2, ] * (1 + 5e-13)
  cleaned <- validate_StochMat(roundoff)
  chain <- MarkovChain(rep(1 / 3, 3), P)
  chain$P <- roundoff
  for (covariance in list(vcov_p, function(x) vcov_pk(x, 3), vcov_stat,
                          function(x) vcov_pk_sl(x, 3))) {
    expect_equal(covariance(roundoff), covariance(cleaned), tolerance = 1e-14)
    expect_equal(covariance(chain), covariance(cleaned), tolerance = 1e-14)
    expect_true(all(is.finite(covariance(roundoff))))
  }
})
