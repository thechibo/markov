# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate Markov Chain Trajectories                                       ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Simulate Markov Chain Trajectories
#'
#' @description Generate one or more independent sample paths from a
#' finite-state, time-homogeneous Markov chain.
#'
#' @param m A single positive integer giving the number of observations per
#'   trajectory, including the initial state.
#' @param x A `MarkovChain` object.
#' @param n A single positive integer giving the number of trajectories.
#'
#' @details
#' Each trajectory starts with a fresh draw from the initial distribution
#' `x$a`. Each subsequent state is drawn from the row of `x$P` for the
#' current state. Simulating multiple paths with a fixed random seed gives
#' the same draws as successive single-path calls with that seed.
#'
#' @return
#' If `n = 1`, an ordered factor of length `m` with one level per state in
#' chain order. If `n > 1`, an `n` by `m` character matrix of state labels,
#' with one trajectory per row and observations in time order across columns.
#' Pass a matrix with at least two columns directly to [fit()] to pool its
#' within-trajectory transitions.
#'
#' @seealso [MarkovChain()], [fit()], [stat()]
#'
#' @export
#'
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
#' x <- MarkovChain(a = c(0.5, 0.5), P = P)
#' set.seed(1)
#' sim(m = 5, x = x)
#' paths <- sim(m = 20, x = x, n = 10)
#' fit(paths)
sim <- function(m, x, n = 1) {
  if (length(m) != 1L || !is.null(dim(m)) || !is_positive_integer(m)) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }
  if (length(n) != 1L || !is.null(dim(n)) || !is_positive_integer(n)) {
    stop("`n` must be a single positive integer.", call. = FALSE)
  }
  m <- as.integer(m)
  n <- as.integer(n)
  x <- validate_MarkovChain(x)
  labels <- get_state_names(x)
  states <- nstates(x)

  draw_path <- function() {
    idx <- integer(m)
    idx[1L] <- sample.int(states, size = 1L, prob = x$a)
    if (m >= 2L) {
      for (time in seq.int(2L, m)) {
        idx[time] <- sample.int(states, size = 1L, prob = x$P[idx[time - 1L], ])
      }
    }
    labels[idx]
  }

  if (n == 1L) {
    return(ordered(draw_path(), levels = labels))
  }

  paths <- vapply(seq_len(n), function(i) draw_path(), character(m))
  matrix(paths, nrow = n, ncol = m, byrow = TRUE)
}
