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
#' @param up Nonempty integer-valued numeric vector of state positions, or
#'   character vector of state names, giving the operational states. Entries
#'   must be distinct and select a proper subset of the states. Input order
#'   is ignored; states are used in their order in `x`.
#'
#' @details
#' The operational subchain must be transient, so that
#' \eqn{I_r-P_{\mathcal{U},\mathcal{U}}} is invertible. The asymptotic
#' variance additionally requires an irreducible chain and treats the initial
#' distribution as fixed.
#'
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
#' @inherit var_mttf references
#'
#' @seealso [avail()] for the time-dependent availability function.
#'
#' @export
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
#' x <- MarkovChain(c(1, 0), P)
#' mttf(x, up = 1)
#' var_mttf(x, up = 1)
mttf <- function(x, up) {

  x <- validate_MarkovChain(x)
  up <- validate_mttf_up(up, get_state_names(x))

  r <- length(up)
  vec_ones <- rep(1, r)
  B <- solve(I_minus_P(x$P, up))
  as.numeric(x$a[up] %*% B %*% vec_ones)

}

#' Asymptotic Variance of the Mean-Time-to-Failure Estimator
#'
#' @description
#' Computes the asymptotic variance of the plug-in estimator of the mean time
#' to failure (MTTF) for a finite-state Markov chain.
#'
#' @param x An object of class \code{MarkovChain}.
#' @param up Nonempty integer-valued numeric vector of state positions, or
#'   character vector of state names, giving the operational states. Entries
#'   must be distinct and select a proper subset of the states. Input order
#'   is ignored; states are used in their order in `x`.
#'
#' @details
#' The operational subchain must be transient, so that
#' \eqn{I_r-P_{\mathcal{U},\mathcal{U}}} is invertible. The asymptotic
#' variance additionally requires an irreducible chain and treats the initial
#' distribution as fixed.
#'
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
#' Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026).
#' \emph{Matrix asymptotic calculus for plug-in maximum likelihood estimators in finite Markov chains}.
#' arXiv preprint arXiv:2607.17161. \url{https://arxiv.org/abs/2607.17161}.
#' See Proposition 11 for the MTTF differential and asymptotic variance.
#'
#' Sadek, A., and Limnios, N. (2002).
#' \emph{Asymptotic properties for maximum likelihood estimators for reliability and failure rates of Markov chains}.
#' Communications in Statistics - Theory and Methods, 31(10), 1837--1861.
#' \doi{10.1081/STA-120014916}.
#'
#' @seealso [mttf()]
#'
#' @export
var_mttf <- function(x, up) {

  x <- validate_MarkovChain(x)
  up <- validate_mttf_up(up, get_state_names(x))

  s <- nstates(x)
  r <- length(up)
  vec_ones <- rep(1, r)
  state_names <- get_state_names(x)
  Nu <- solve(I_minus_P(x$P, up))

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
#' `var_mttf_sl()` returns a numeric scalar containing the asymptotic variance.
#' `deriv_psi()` returns a numeric vector of length `length(up)^2`, containing
#' the derivatives with respect to the operational transition block in
#' row-wise order, with operational states ordered as in `x`.
#'
#' @references
#' Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026).
#' \emph{Matrix asymptotic calculus for plug-in maximum likelihood estimators in finite Markov chains}.
#' arXiv preprint arXiv:2607.17161. \url{https://arxiv.org/abs/2607.17161}.
#' See Proposition 11 for the MTTF differential and asymptotic variance.
#'
#' Sadek, A., and Limnios, N. (2002).
#' \emph{Asymptotic properties for maximum likelihood estimators for reliability and failure rates of Markov chains}.
#' Communications in Statistics - Theory and Methods, 31(10), 1837--1861.
#' \doi{10.1081/STA-120014916}.
#'
#' @seealso [mttf()], [var_mttf()]
#'
#' @keywords internal
var_mttf_sl <- function(x, up) {

  x <- validate_MarkovChain(x)
  up <- validate_mttf_up(up, get_state_names(x))

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

  x <- validate_MarkovChain(x)
  up <- validate_mttf_up(up, get_state_names(x))

  s1 <- length(up)
  psi_prime <- numeric(s1 ^ 2)
  col_ones <- matrix(1, ncol = 1, nrow = s1)
  Lambda <- solve(I_minus_P(x$P, up))

  for (i in 1:(s1 ^ 2)) {
    mat <- matrix(cron(i, s1 ^ 2), nrow = s1, ncol = s1, byrow = TRUE)
    psi_prime[i] <- x$a[up] %*% Lambda %*% mat %*% Lambda %*% col_ones
  }

  psi_prime

}

# Normalize positions and names to the chain's order so that the derivative
# and covariance mask use the same coordinates.
validate_mttf_up <- function(up, state_names) {
  if (!is.null(dim(up)) || !(is.numeric(up) || is.character(up)) ||
      length(up) == 0L || anyNA(up)) {
    stop("`up` must be a nonempty vector of state names or integer indices.",
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
  if (anyDuplicated(up)) {
    stop("`up` must contain distinct states.", call. = FALSE)
  }
  if (length(up) == length(state_names)) {
    stop("`up` must be a proper subset of the states.", call. = FALSE)
  }
  state_names[state_names %in% up]
}
