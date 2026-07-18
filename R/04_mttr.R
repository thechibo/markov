# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mean Time to Recovery                                                     ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Mean Time to Recovery
#'
#' @description
#' Computes the mean time to recovery of a discrete-time, finite-state Markov
#' chain. Revocery is defined as the first time the chain leaves the set of
#' operational states specified by `down`.
#'
#' @param x An object of class \code{MarkovChain}.
#' @param down Integer or character vector giving the down (non-operational)
#' states.
#'
#' @details
#' Let \eqn{\mathcal{D}} denote the set of non-operational states,
#' \eqn{P_{\mathcal{D},\mathcal{D}}} be the submatrix of transition
#' probabilities between non-operational states, \eqn{d=|\mathcal{D}|}, and
#' \eqn{N=(I_d-P_{\mathcal{D},\mathcal{D}})^{-1}}. The mean time to recovery is
#' computed as
#' \deqn{
#'   \psi
#'   =
#'   a_{\mathcal{D}}^\mathsf{T}
#'   N
#'   \mathbf{1}_d,
#' }
#' where \eqn{a_{\mathcal{D}}} is the restriction of the initial distribution
#' to the non-operational states and \eqn{\mathbf{1}_d} is a vector of ones.
#'
#' Let \eqn{\Sigma_{\mathcal{D}}} be the submatrix of \eqn{\Sigma_{\mathbf p}}
#' corresponding to transitions between non-operational states. Then:
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
#'   \Sigma_{\mathcal{D}}
#'   D_\psi^{\mathsf T},
#' }
#' and
#' \deqn{
#'   D_\psi
#'   =
#'   \left(
#'     a_{\mathcal{D}}^{\mathsf T}N
#'   \right)
#'   \mathbin{\otimes}
#'   \left(
#'     \mathbf{1}_d^{\mathsf T}N^{\mathsf T}
#'   \right).
#' }
#'
#' The MTTR \eqn{\psi} is returned by this function.
#' The asymptotic variance \eqn{\sigma_\psi^2} is returned by [var_mttr()].
#'
#' @return
#' A numeric scalar giving the expected number of discrete-time transitions
#' before the chain first enters a state outside `down`.
#'
#' @seealso [mttf()], [var_mttr()]
#'
#' @export
#' @inherit MarkovChain examples
mttr <- function(x, down) {
  mttf(x, down)
}

#' Asymptotic Variance of the Mean-Time-to-Recovery Estimator
#'
#' @description
#' Computes the asymptotic variance of the plug-in estimator of the mean time
#' to recovery (MTTR) for a finite-state Markov chain.
#'
#' @param x An object of class \code{MarkovChain}.
#' @param down Integer or character vector giving the down (non-operational)
#' states.
#'
#' @details
#' Let \eqn{\mathcal{D}} denote the set of non-operational states,
#' \eqn{P_{\mathcal{D},\mathcal{D}}} be the submatrix of transition
#' probabilities between non-operational states, \eqn{d=|\mathcal{D}|}, and
#' \eqn{N=(I_d-P_{\mathcal{D},\mathcal{D}})^{-1}}. The mean time to recovery is
#' computed as
#' \deqn{
#'   \psi
#'   =
#'   a_{\mathcal{D}}^\mathsf{T}
#'   N
#'   \mathbf{1}_d,
#' }
#' where \eqn{a_{\mathcal{D}}} is the restriction of the initial distribution
#' to the non-operational states and \eqn{\mathbf{1}_d} is a vector of ones.
#'
#' Let \eqn{\Sigma_{\mathcal{D}}} be the submatrix of \eqn{\Sigma_{\mathbf p}}
#' corresponding to transitions between non-operational states. Then:
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
#'   \Sigma_{\mathcal{D}}
#'   D_\psi^{\mathsf T},
#' }
#' and
#' \deqn{
#'   D_\psi
#'   =
#'   \left(
#'     a_{\mathcal{D}}^{\mathsf T}N
#'   \right)
#'   \mathbin{\otimes}
#'   \left(
#'     \mathbf{1}_d^{\mathsf T}N^{\mathsf T}
#'   \right).
#' }
#'
#' The MTTR \eqn{\psi} is returned by [mttr()].
#' The asymptotic variance \eqn{\sigma_\psi^2} is returned by this function.
#'
#' @return
#' A numeric scalar giving the expected number of discrete-time transitions
#' before the chain first enters a state outside `down`.
#'
#' @seealso [mttf()], [var_mttr()]
#'
#' @export
#' @inherit MarkovChain examples
var_mttr <- function(x, down) {
  var_mttf(x, down)
}
