# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate a Markov Chain Trajectory                                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Simulate a Markov Chain Trajectory
#'
#' @description Generate a single sample path of length `m` from a finite-state,
#' time-homogeneous Markov chain.
#'
#' @param m Positive integer. The desired trajectory length.
#' @param x A `MarkovChain` object.
#'
#' @details
#' The initial state is drawn from the chain's initial distribution `a`, and
#' subsequent states are generated iteratively using the transition matrix `P`.
#'
#' @return
#' A factor of length `m` giving the visited states in order. Levels correspond
#' to the state labels (from `names(x$a)` or `rownames(x$P)` if available).
#'
#' @seealso [MarkovChain()], [mc_estim()], [mc_stat()]
#'
#' @export
#'
#' @inherit MarkovChain examples
mc_sim <- function(m, x) {

  # Preliminaries
  m <- is_positive_integer(m)
  x <- validate_MarkovChain(x)
  st_names <- get_state_names(x)
  idx <- integer(m)
  states <- nstates(x)

  # simulate the initial state
  idx[1] <- sample(states, size = 1, prob = x$a)

  # simulate the rest
  if (m >= 2L) {
    for (i in seq.int(2L, m)) {
      idx[i] <- sample(states, size = 1L, prob = x$P[idx[i - 1L], ])
    }
  }

  factor(st_names[idx], levels = st_names, ordered = TRUE)
}
