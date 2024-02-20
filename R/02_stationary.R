# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stationary Distribution
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Stationary Distribution
#'
#' @param x an object of class `MarkovChain`.
#'
#' @return A stochastic vector with the chain's stationary distribution.
#' @export
#'
#' @inherit MarkovChain examples
stationary <- function(x) {

  r <- eigen(t(x$P))
  rvec <- r$vectors[, which.min(abs(r$values - 1))]
  rvec / sum(rvec)

  StochVec(rvec / sum(rvec))

}
