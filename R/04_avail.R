# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Availability                                                              ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Availability of a Markov Chain
#'
#' @description
#' Compute the discrete-time availability function of a system whose
#' degradation/repair dynamics are modeled by a finite-state Markov chain.
#' Availability at time \eqn{t} is defined as
#' \deqn{
#'   A(t) = \Pr\{ X_t \in \mathcal{U} \}
#' }
#' where \eqn{X_t} is the state of the chain at time \eqn{t} and
#' \eqn{\mathcal{U}} is the set of \emph{up} (operational) states.
#'
#' @param x An object of class \code{MarkovChain}.
#' @param time A positive integer vector giving the discrete time points
#'   \eqn{t} at which availability \eqn{A(t)} is evaluated.
#' @param up Integer or character vector giving the up (operational) states.
#'
#' @details
#' For each time \eqn{t}, the function computes
#' \deqn{
#'   A(t) = a P^t \mathbf{1}_{\mathcal{U}},
#' }
#' where \eqn{a} is the initial distribution, \eqn{P} is the transition
#' matrix, and \eqn{\mathbf{1}_{\mathcal{U}}} is the indicator vector of the
#' operational states.
#'
#' Duplicate values in `up` are ignored.
#'
#' @return
#' A numeric vector of the same length as `time`, where the \eqn{j}-th element
#' is the probability that the system is available (in an up state) at
#' `time[j]`.
#'
#' @seealso [MarkovChain()], [Ppower()]
#'
#' @export
#' @inherit MarkovChain examples
avail <- function(x, time, up) {

  # Preliminary Checks
  if (!inherits(x, "MarkovChain")) {
    stop("`x` must be a `MarkovChain` object.", call. = FALSE)
  }
  if (!is_positive_integer(time)) {
    stop("`time` must be a positive integer vector.", call. = FALSE)
  }
  up <- as.character(up)
  if (!all(up %in% get_state_names(x))) {
    stop("`up` must be a vector of state names or indices.", call. = FALSE)
  }

  vec_ones <- as.numeric(get_state_names(x) %in% up)
  y <- numeric(length(time))
  for (i in seq_along(time)) {
    y[i] <- as.numeric(x$a %*% Ppower(x$P, time[i]) %*% vec_ones)
  }
  y
}

#' Variance of availability
#'
#' @description
#' Compute the variance of the system availability at discrete time \eqn{n}
#' for a finite-state ergodic Markov chain.
#'
#' @param x An object of class \code{MarkovChain}.
#' @param time Integer scalar. Discrete time step \eqn{time} at which the
#' variance of availability is evaluated.
#' @param up Integer or character vector giving the up (operational) states.
#'
#' @return
#' A numeric scalar: the variance of the availability \eqn{A(n)} at time \eqn{n}
#' as defined in Sadek and Limnios (2002).
#'
#' @references
#' A. Sadek and N. Limnios (2002).
#' \emph{Asymptotic properties of estimators for nonhomogeneous Markov chains}.
#' Statistics & Probability Letters, 60(4), 1839–1850.
#' See Theorem 5.
#'
#' @seealso [avail()], [vcov_p()], [vcov_pk()]
#'
#' @export
vcov_avail <- function(x, time, up) {

  # Preliminary Checks
  if (!inherits(x, "MarkovChain")) {
    stop("`x` must be a `MarkovChain` object.", call. = FALSE)
  }
  if (!is_positive_integer(time)) {
    stop("`time` must be a positive integer vector.", call. = FALSE)
  }
  up <- as.character(up)
  if (!all(up %in% get_state_names(x))) {
    stop("`up` must be a vector of state names or indices.", call. = FALSE)
  }

  s <- matrix(0, ncol = nstates(x) ^ 2, nrow = 1)
  vec_ones <- matrix(as.numeric(get_state_names(x) %in% up), nrow = 1)

  for (i in 1:time) {
    tmp <- fastmatrix::kronecker.prod(x$a %*% Ppower(x$P, i - 1),
                                      t(Ppower(x$P, time - i) %*% t(vec_ones)))
    s <- s + tmp
  }

  as.numeric(s %*% vcov_p(x) %*% t(s))

}

#' @keywords internal
vcov_avail_sl <- function(x, time, up) {

  S <- vcov_pk(x, time)
  k <- nstates(x)
  down <- setdiff(1:k, up)

  for (i in seq_along(x$a)) {
    m <- (i - 1) * k + down
    S[, m] <- 0
  }

  init <- rep(c(x$a[up], rep(0, length(down))), each = length(x$a))
  as.numeric(t(init) %*% S %*% init)

}
