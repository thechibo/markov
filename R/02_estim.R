# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Estimate Markov Chain
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Markov Chain Parameter Estimation
#'
#' @param x vector. A trajectory of a Markov Chain.
#'
#' @return An object of class `MarkovChain`.
#' @export
#'
#' @inherit MarkovChain examples
estim <- function(x) {

  # Preliminaries
  states <- unique(x)
  if (is.numeric(states) || is.ordered(states)) {
    states <- sort(states)
  }
  s <- length(states)

  P <- matrix(0, nrow = s, ncol = s, byrow = TRUE)
  dimnames(P) <- list(states, states)

  # Estimate the transition probabilities
  for (t in 1:(length(x) - 1)) {
    P[x[t], x[t + 1]] <- P[x[t], x[t + 1]] + 1
  }

  it <- 0
  while (any(rowSums(P) != 1) & it < 10) {
    P <- P / rowSums(P)
    it <- it + 1
  }
  if (it > 2) {
    warning("The transition matrix was rescaled due to numerical errors.")
  }

  # Initial state
  # a <- rep(0, length = s)
  # a[x[1]] <- 1
  a <- rep(1/s, length = s)

  MarkovChain(a = StochVec(a), P = StochMat(P))

}
