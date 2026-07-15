# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Availability                                                              ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Availability of a Markov chain
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
#' Equivalently, for each time \eqn{t} in `time`, the function computes
#' \deqn{
#'   A(t) = a P^t \mathbf{1}_{\mathcal{U}},
#' }
#' where \eqn{P = x$P} is the transition matrix and
#' \eqn{\mathbf{1}_{\mathcal{U}}} is the indicator vector of up states.
#'
#' @param x An object of class \code{MarkovChain}.
#' @param time Integer vector. Discrete time steps \eqn{t} at which to evaluate
#'   the availability \eqn{A(t)}.
#'
#' @return
#' A numeric vector of the same length as `time`, where the \eqn{j}-th element
#' is the probability that the system is available (in an up state) at time
#' `time[j]`.
#'
#' @seealso [MarkovChain()], [mat_gamma()]
#'
#' @export
#' @inherit MarkovChain examples
avail <- function(x, time, up) {

  k <- nstates(x)
  nup <- length(up)
  vec_ones <- as.numeric(1:k %in% up)

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
#' for a finite-state ergodic Markov chain, using the asymptotic covariance
#' structure described in Theorem 5 of Sadek and Limnios (2002).
#'
#' @param x An object of class \code{MarkovChain}.
#' @param time Integer scalar. Discrete time step \eqn{time} at which the
#' variance of availability is evaluated.
#' @param up Integer vector giving the indices of up (operational) states.
#'
#' @return
#' A numeric scalar: the variance of the availability \eqn{A(n)} at time \eqn{n}
#' as defined in Sadek and Limnios (2002).
#'
#' @references
#' A. Sadek and N. Limnios (2002).
#' \emph{Asymptotic properties of estimators for nonhomogeneous Markov chains}.
#' Statistics & Probability Letters, 60(4), 1839–1850.
#'
#' @seealso [avail()], [mat_gamma()], [mat_sigma()]
#'
#' @export
var_avail <- function(x, time, up) {

  k <- nstates(x)
  s <- matrix(0, ncol = k ^ 2, nrow = 1)

  vec_ones <- matrix(0, ncol = k, nrow = 1)
  vec_ones[up] <- 1

  for (i in 1:time) {
    tmp <- fastmatrix::kronecker.prod(x$a %*% Ppower(x$P, i - 1),
                                      t(Ppower(x$P, time - i) %*% t(vec_ones)))
    s <- s + tmp
  }

  as.numeric(s %*% mat_gamma(x) %*% t(s))

}

#' @keywords internal
var_avail_sl <- function(x, time, up) {

  S <- mat_sigma(x, time)
  k <- nstates(x)
  down <- setdiff(1:k, up)

  for (i in seq_along(x$a)) {
    m <- (i - 1) * k + down
    S[, m] <- 0
  }

  init <- rep(c(x$a[up], rep(0, length(down))), each = length(x$a))
  as.numeric(t(init) %*% S %*% init)

}
