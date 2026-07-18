# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mean Time to Failure                                                      ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Mean Time to Failure
#'
#' @description
#' Computes the mean time to failure of a discrete-time, finite-state Markov
#' chain. Failure is defined as the first time the chain leaves the set of
#' operational states specified by `up`.
#'
#' @param x An object of class \code{MarkovChain}.
#' @param up Integer or character vector giving the up (operational) states.
#'
#' @details
#' Let \eqn{\mathcal{U}} denote the set of operational states,
#' \eqn{P_{\mathcal{U},\mathcal{U}}} be the submatrix of transition
#' probabilities between operational states, \eqn{r=|\mathcal{U}|}, and
#' \eqn{N=(I_r-P_{\mathcal{U},\mathcal{U}})^{-1}}. The mean time to failure is
#' computed as
#' \deqn{
#'   \psi
#'   =
#'   a_{\mathcal{U}}^\mathsf{T}
#'   N
#'   \mathbf{1}_r,
#' }
#' where \eqn{a_{\mathcal{U}}} is the restriction of the initial distribution
#' to the operational states and \eqn{\mathbf{1}_r} is a vector of ones.
#'
#' Let \eqn{\Sigma_{\mathcal{U}}} be the submatrix of \eqn{\Sigma_{\mathbf p}}
#' corresponding to transitions between operational states. Then:
#' \deqn{
#'   \sqrt{n}
#'   \left\{
#'     \widehat{\psi}_n-\psi
#'   \right\}
#'   \mathrel{\mathop{\longrightarrow}\limits^{d}}
#'   \mathcal{N}
#'   \left(
#'     0,
#'     \sigma_\psi^2
#'   \right),
#' }
#' where
#' \deqn{
#'   \sigma_\psi^2
#'   =
#'   D_\psi
#'   \Sigma_{\mathcal{U}}
#'   D_\psi^{\mathsf T},
#' }
#' and
#' \deqn{
#'   D_\psi
#'   =
#'   \left(
#'     a_{\mathcal{U}}^{\mathsf T}N
#'   \right)
#'   \mathbin{\otimes}
#'   \left(
#'     \mathbf{1}_r^{\mathsf T}N^{\mathsf T}
#'   \right).
#' }
#'
#' The MTTF \eqn{\psi} is returned by this function.
#' The asymptotic variance \eqn{\sigma_\psi^2} is returned by [var_mttf()].
#'
#' @return
#' A numeric scalar giving the expected number of discrete-time transitions
#' before the chain first enters a state outside `up`.
#'
#' @seealso [avail()] for the time-dependent availability function.
#'
#' @export
#' @inherit MarkovChain examples
mttf <- function(x, up) {

  # Preliminary Checks
  if (!inherits(x, "MarkovChain")) {
    stop("`x` must be a `MarkovChain` object.", call. = FALSE)
  }
  up <- as.character(up)
  if (!all(up %in% get_state_names(x))) {
    stop("`up` must be a vector of state names or indices.", call. = FALSE)
  }

  r <- length(up)
  vec_ones <- rep(1, r)
  B <- solve(diag(r) - x$P[up, up])
  as.numeric(x$a[up] %*% B %*% vec_ones)

}

#' Asymptotic Variance of the Mean-Time-to-Failure Estimator
#'
#' @description
#' Computes the asymptotic variance of the plug-in estimator of the mean time
#' to failure (MTTF) for a finite-state Markov chain.
#'
#' @param x An object of class \code{MarkovChain}.
#' @param up Integer or character vector giving the up (operational) states.
#'
#' @details
#' Let \eqn{\mathcal{U}} denote the set of operational states,
#' \eqn{P_{\mathcal{U},\mathcal{U}}} be the submatrix of transition
#' probabilities between operational states, \eqn{r=|\mathcal{U}|}, and
#' \eqn{N=(I_r-P_{\mathcal{U},\mathcal{U}})^{-1}}. The mean time to failure is
#' computed as
#' \deqn{
#'   \psi
#'   =
#'   a_{\mathcal{U}}^\mathsf{T}
#'   N
#'   \mathbf{1}_r,
#' }
#' where \eqn{a_{\mathcal{U}}} is the restriction of the initial distribution
#' to the operational states and \eqn{\mathbf{1}_r} is a vector of ones.
#'
#' Let \eqn{\Sigma_{\mathcal{U}}} be the submatrix of \eqn{\Sigma_{\mathbf p}}
#' corresponding to transitions between operational states. Then:
#' \deqn{
#'   \sqrt{n}
#'   \left\{
#'     \widehat{\psi}_n-\psi
#'   \right\}
#'   \mathrel{\mathop{\longrightarrow}\limits^{d}}
#'   \mathcal{N}
#'   \left(
#'     0,
#'     \sigma_\psi^2
#'   \right),
#' }
#' where
#' \deqn{
#'   \sigma_\psi^2
#'   =
#'   D_\psi
#'   \Sigma_{\mathcal{U}}
#'   D_\psi^{\mathsf T},
#' }
#' and
#' \deqn{
#'   D_\psi
#'   =
#'   \left(
#'     a_{\mathcal{U}}^{\mathsf T}N
#'   \right)
#'   \mathbin{\otimes}
#'   \left(
#'     \mathbf{1}_r^{\mathsf T}N^{\mathsf T}
#'   \right).
#' }
#'
#' The MTTF \eqn{\psi} is returned by [mttf()].
#' The asymptotic variance \eqn{\sigma_\psi^2} is returned by this function.
#'
#' @return
#' A numeric scalar containing the asymptotic variance \eqn{\sigma_\psi^2} of
#' the MTTF estimator.
#'
#' @references
#' GSI
#'
#' A. Sadek and N. Limnios (2002).
#' \emph{Asymptotic properties of estimators for nonhomogeneous Markov chains}.
#' Statistics & Probability Letters, 60(4), 1839–1850.
#' See Theorem 8 (p. 1852) for the MTTF.
#'
#' @seealso [mttf()]
#'
#' @export
var_mttf <- function(x, up) {

  # Preliminary Checks
  if (!inherits(x, "MarkovChain")) {
    stop("`x` must be a `MarkovChain` object.", call. = FALSE)
  }
  up <- as.character(up)
  if (!all(up %in% get_state_names(x))) {
    stop("`up` must be a vector of state names or indices.", call. = FALSE)
  }

  s <- nstates(x)
  r <- length(up)
  vec_ones <- rep(1, r)
  state_names <- get_state_names(x)
  Nu <- solve(diag(r) - x$P[up, up])

  M <- kp(x$a[up] %*% Nu, t(vec_ones) %*% t(Nu))

  upmat <- matrix(FALSE, nrow = s, ncol = s,
                  dimnames = list(state_names, state_names))
  upmat[up, up] <- TRUE
  upvec <- as.vector(upmat)

  as.vector(M %*% vcov_p(x)[upvec, upvec] %*% t(M))

}

#' Reference Implementation of the Asymptotic Variance of MTTF
#'
#' @description
#' Computes the asymptotic variance of the mean-time-to-failure estimator
#' using a component-wise derivative calculation.
#'
#' @inheritParams var_mttf
#'
#' @details
#' Unlike [var_mttf()], which uses a matrix and Kronecker-product
#' implementation, this function follows the corresponding scalar or
#' component-wise expressions. It is used for internal testing only.
#'
#' @return
#' A numeric scalar containing the asymptotic variance of the MTTF estimator.
#'
#' @references
#' A. Sadek and N. Limnios (2002).
#' \emph{Asymptotic properties of estimators for nonhomogeneous Markov chains}.
#' Statistics & Probability Letters, 60(4), 1839–1850.
#' See Theorem 8 (p. 1852) for the MTTF.
#'
#' @seealso [mttf()], [var_mttf()]
#'
#' @keywords internal
var_mttf_sl <- function(x, up) {

  psi_prime <- deriv_psi(x, up)
  s <- nstates(x)
  state_names <- get_state_names(x)

  upmat <- matrix(FALSE, nrow = s, ncol = s,
                  dimnames = list(state_names, state_names))
  upmat[up, up] <- TRUE
  upvec <- as.vector(upmat)

  as.vector(t(psi_prime) %*% vcov_p(x)[upvec, upvec] %*% psi_prime)

}

#' @rdname var_mttf_sl
deriv_psi <- function(x, up) {

  # Preliminary Checks
  if (!inherits(x, "MarkovChain")) {
    stop("`x` must be a `MarkovChain` object.", call. = FALSE)
  }
  up <- as.character(up)
  if (!all(up %in% get_state_names(x))) {
    stop("`up` must be a vector of state names or indices.", call. = FALSE)
  }

  s1 <- length(up)
  psi_prime <- numeric(s1 ^ 2)
  col_ones <- matrix(1, ncol = 1, nrow = s1)
  a_up <- x$a[up]
  Lambda <- solve(diag(s1) - x$P[up, up]) # inverse of (I - P11)

  for (i in 1:(s1 ^ 2)) {
    mat <- matrix(cron(i, s1 ^ 2), nrow = s1, ncol = s1, byrow = TRUE)
    psi_prime[i] <- x$a[up] %*% Lambda %*% mat %*% Lambda %*% col_ones
  }

  psi_prime

}
