# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Estimate Markov Chain Parameters                                          ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Estimate a Markov chain from an observed state sequence
#'
#' @description
#' Estimates a time-homogeneous, discrete-time Markov chain from an observed
#' sequence of states.
#'
#' @param x A vector of observed states in time order. Can be numeric, character,
#'   factor, or ordered. Must have length at least 2 and contain no `NA`s.
#' @param init Initial distribution to attach to the fitted Markov chain.
#'   Must be exactly one of:
#'   \describe{
#'     \item{`"uniform"`}{Uniform distribution over the observed state space.}
#'     \item{`"realized"`}{Degenerate distribution concentrated on the first observed state `x[1]`. This is the default choice.}
#'     \item{`"stationary"`}{Stationary distribution implied by the estimated transition matrix (via `stat()`).}
#'   }
#'
#' @details
#' The state space is defined as `unique(x)`. If `x` is numeric or ordered, the
#' states are sorted; otherwise, the first-appearance order is used.
#'
#' The function computes the empirical transition frequencies and normalizes by
#' row to obtain a transition probability matrix.
#'
#' If a state never appears as an origin (i.e., it is never observed at times
#' `t = 1, ..., n-1`), the corresponding transition row has no empirical counts.
#' In that case the function treats that state as absorbing (probability 1 on the
#' diagonal) and emits a warning.
#'
#' A numerical cleanup step will rescale rows if their sums deviate from 1 by
#' by more than 1e-12 (tolerance), emitting a warning.
#'
#' @return A `MarkovChain` object.
#'
#' @seealso [MarkovChain()], [sim()], [stat()]
#'
#' @inherit MarkovChain examples
#' @export
fit <- function(x, init = c("realized", "uniform", "stationary")) {

  # init must be a single choice
  init <- match.arg(init)

  # Basic checks
  if (length(x) < 2) {
    stop("`x` must have length >= 2 to estimate transitions.")
  }
  if (anyNA(x)) {
    stop("`x` contains NA values.
         Please remove/impute them before calling `fit()`.")
  }

  # States (keep stable ordering for numeric/ordered)
  states <- unique(x)
  if (is.numeric(states) || is.ordered(states)) {
    states <- sort(states)
  }
  s <- length(states)

  # Map x to factor levels to ensure consistent indexing
  xf <- factor(x, levels = states)

  # Transition counts
  from <- xf[-length(xf)]
  to   <- xf[-1L]

  P <- matrix(0, nrow = s, ncol = s, dimnames = list(states, states))
  counts <- table(from, to)
  P[rownames(counts), colnames(counts)] <- as.numeric(counts)

  # Row-normalize to transition probabilities
  rs <- rowSums(P)

  # Handle states that never appear as an origin (row sum = 0)
  if (any(rs == 0)) {
    zero_states <- states[rs == 0]
    warning(
      "Some states never appear as an origin in `x`;
      treating them as absorbing states: ",
      paste(zero_states, collapse = ", ")
    )
    for (st in zero_states) P[st, st] <- 1
    rs <- rowSums(P)
  }

  P <- sweep(P, 1, rs, "/")

  # Numerical cleanup
  P[abs(P) < .Machine$double.eps] <- 0
  rs2 <- rowSums(P)
  if (any(abs(rs2 - 1) > 1e-12)) {
    P <- sweep(P, 1, rs2, "/")
    warning("The transition matrix was rescaled due to numerical errors.")
  }
  P <- StochMat(P)

  # Initial distribution
  if (init == "uniform") {
    a <- rep(1 / s, s)
    names(a) <- states
  } else if (init == "realized") {
    a <- setNames(rep(0, s), states)
    a[as.character(xf[1L])] <- 1
  } else { # init == "stationary"
    a <- stat(P)
    a <- setNames(as.numeric(a), states)
  }

  MarkovChain(a = StochVec(a), P = P)
}
