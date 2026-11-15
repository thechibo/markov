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
#' @param time A nonempty positive integer vector giving the discrete time points
#'   \eqn{t} at which availability \eqn{A(t)} is evaluated.
#' @param up Integer-valued numeric vector of state positions or character
#'   vector of state names giving the operational states. Missing values are
#'   not allowed. Duplicates are ignored; empty and full state sets are valid.
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
#' Powers are computed by repeated real matrix multiplication with [Ppowers()].
#'
#' @return
#' A numeric vector of the same length as `time`, where the \eqn{j}-th element
#' is the probability that the system is available (in an up state) at
#' `time[j]`.
#'
#' @seealso [MarkovChain()], [Ppowers()], [vcov_avail()]
#'
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
#' x <- MarkovChain(c(0.4, 0.6), P)
#' avail(x, time = 1:3, up = 1L)
#' vcov_avail(x, time = 3L, up = 1L)
#'
#' @export
avail <- function(x, time, up) {

  x <- validate_MarkovChain(x)
  if (length(time) == 0L || !is.null(dim(time)) ||
      !is_positive_integer(time)) {
    stop("`time` must be a nonempty positive integer vector.", call. = FALSE)
  }
  up <- validate_avail_up(up, get_state_names(x))

  vec_ones <- as.numeric(get_state_names(x) %in% up)
  powers <- Ppowers(x$P, max(time))
  y <- numeric(length(time))
  for (i in seq_along(time)) {
    y[i] <- as.numeric(x$a %*% powers[[time[i] + 1L]] %*% vec_ones)
  }
  y
}

#' Variance of availability
#'
#' @description
#' Compute the pointwise asymptotic variance of the plug-in availability
#' estimator for an irreducible finite-state Markov chain, treating the
#' initial distribution as fixed.
#'
#' @param x An object of class \code{MarkovChain}.
#' @param time A single positive integer giving the time at which the
#' asymptotic variance is evaluated.
#' @param up Integer-valued numeric vector of state positions or character
#'   vector of state names giving the operational states. Missing values are
#'   not allowed. Duplicates are ignored; empty and full state sets are valid.
#'
#' @details
#' This function evaluates the diagonal entry of Equation (29) for the
#' requested time, using the operational-state indicator specified in
#' Proposition 9 of Gavrilopoulos, Trevezas, and Votsi (2026). The derivative
#' is assembled as a sum of Kronecker products, using real transition powers
#' computed once by [Ppowers()]. No eigendecomposition or inverse eigenvector
#' matrix is required.
#'
#' @return
#' A numeric scalar containing the variance in the \eqn{\sqrt{m}} limit
#' of the availability estimator, where \eqn{m} is the sample size.
#' Divide by \eqn{m} to approximate the variance of the estimator itself.
#' This is a pointwise variance, not a joint covariance across times.
#'
#' @references
#' Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026).
#' \emph{Matrix asymptotic calculus for plug-in maximum likelihood estimators in finite Markov chains}.
#' arXiv preprint arXiv:2607.17161. \url{https://arxiv.org/abs/2607.17161}.
#'
#' @seealso [avail()], [vcov_p()], [vcov_pk()]
#'
#' @inherit avail examples
#' @export
vcov_avail <- function(x, time, up) {

  x <- validate_MarkovChain(x)
  validate_avail_time(time)
  up <- validate_avail_up(up, get_state_names(x))

  s <- matrix(0, ncol = nstates(x) ^ 2, nrow = 1)
  vec_ones <- matrix(as.numeric(get_state_names(x) %in% up), nrow = 1)
  # Ppowers() requires a positive highest exponent; time = 1 only uses P^0.
  powers <- Ppowers(x$P, max(1L, time - 1L))

  for (i in seq_len(time)) {
    tmp <- fastmatrix::kronecker.prod(x$a %*% powers[[i]],
                                      t(powers[[time - i + 1L]] %*% t(vec_ones)))
    s <- s + tmp
  }

  as.numeric(s %*% vcov_p(x) %*% t(s))

}

#' Component-wise Reference for Availability Variance
#'
#' @description Computes the pointwise asymptotic variance using the
#' component-wise covariance expansion of Sadek and Limnios (2002).
#' @inheritParams vcov_avail
#' @details Theorem 5 (p. 1850) takes the variance of the sum of
#' \eqn{a_i f_{ij}^{t}} over all initial states \eqn{i} and operational
#' destination states \eqn{j}. The implementation adds the weighted diagonal
#' variances and twice every unordered pairwise covariance. The covariance
#' entries are computed by [vcov_pk_sl()] using Theorem 4 (p. 1847).
#' Unlike [vcov_avail()], this reference does not construct a Kronecker
#' derivative. Arbitrary operational-state labels replace the paper's
#' convention that the first \eqn{r} states are operational.
#' @return A numeric scalar containing the pointwise asymptotic variance.
#' @references
#' Sadek, A., and Limnios, N. (2002).
#' \emph{Asymptotic properties for maximum likelihood estimators for reliability and failure rates of Markov chains}.
#' Communications in Statistics - Theory and Methods, 31(10), 1837--1861.
#' \doi{10.1081/STA-120014916}.
#' See Theorem 5 (p. 1850) for availability and Theorem 4 (p. 1847)
#' for the covariance of transition powers.
#' @keywords internal
vcov_avail_sl <- function(x, time, up) {
  x <- validate_MarkovChain(x)
  validate_avail_time(time)
  up <- validate_avail_up(up, get_state_names(x))

  S <- vcov_pk_sl(x, time)
  k <- nstates(x)
  destinations <- match(up, get_state_names(x))
  coordinates <- as.vector(t(outer((seq_len(k) - 1L) * k,
                                   destinations, "+")))
  weights <- rep(as.numeric(x$a), each = length(destinations))

  # Theorem 5: diagonal variances plus twice each unordered covariance.
  variance <- sum(weights^2 * diag(S)[coordinates])
  for (j in seq_along(coordinates)) {
    if (j > 1L) {
      previous <- seq_len(j - 1L)
      variance <- variance + 2 * weights[j] * sum(
        weights[previous] * S[coordinates[j], coordinates[previous]]
      )
    }
  }
  as.numeric(variance)
}

validate_avail_time <- function(time) {
  if (length(time) != 1L || !is.null(dim(time)) ||
      !is_positive_integer(time)) {
    stop("`time` must be a single positive integer.", call. = FALSE)
  }
  invisible(NULL)
}

# Availability uses a set: duplicates are ignored, and empty/full sets are valid.
validate_avail_up <- function(up, state_names) {
  if (!is.null(dim(up)) || !(is.numeric(up) || is.character(up)) || anyNA(up)) {
    stop("`up` must be a vector of state names or integer indices.",
         call. = FALSE)
  }
  if (is.numeric(up)) {
    if (any(!is.finite(up)) || any(up != floor(up)) ||
        any(up < 1 | up > length(state_names))) {
      stop("`up` indices must be integers between 1 and the number of states.",
           call. = FALSE)
    }
    up <- state_names[up]
  } else if (!all(up %in% state_names)) {
    stop("`up` contains unknown state names.", call. = FALSE)
  }
  state_names[state_names %in% up]
}

#' Asymptotic Confidence Intervals for Availability
#'
#' Compute two-sided, pointwise normal-approximation confidence intervals using [avail()]
#' and [vcov_avail()].
#' @inheritParams avail
#' @inheritParams ci_rel
#' @inheritSection ci_rel Normal approximation
#' @details
#' The asymptotic normality result is Proposition 9 of Gavrilopoulos, Trevezas, and Votsi (2026).
#' The domain restrictions of [vcov_avail()] apply.
#' @return A data frame with numeric columns `estimate`, `std.error`,
#'   `lower`, and `upper`, giving the point estimate, estimated standard
#'   error, and confidence limits.
#'   There is one row per requested time, preserving order and duplicates.
#' @inherit ci_rel references
#' @seealso [avail()], [vcov_avail()], [ci_rel()]
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE)
#' x <- MarkovChain(c(0.4, 0.6), P)
#' ci_avail(x, time = 1:3, up = 1L, n = 1000)
#' @export
ci_avail <- function(x, time, up, n, level = 0.95) {
  validate_ci_parameters(n, level)
  estimate <- avail(x, time, up)
  times <- unique(time)
  variance <- vapply(times, function(t) vcov_avail(x, t, up),
                     numeric(1), USE.NAMES = FALSE)
  norm_ci(estimate, variance[match(time, times)], n, level)
}
