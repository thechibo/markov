# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mean Time to Recovery                                                     ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Mean Time to Recovery
#'
#' @description
#' Computes the mean time to recovery of a discrete-time, finite-state Markov
#' chain. Recovery is defined as the first time the chain leaves the set of
#' non-operational states specified by `down`.
#'
#' @param x An object of class \code{MarkovChain}.
#' @param down Nonempty integer-valued numeric vector of state positions, or
#'   character vector of state names, giving the non-operational states.
#'   Entries must be distinct and select a proper subset of the states.
#'   Input order is ignored; states are used in their order in `x`.
#'
#' @details
#' These functions use [mttf()] and [var_mttf()] with `down` as the selected
#' state set. The non-operational subchain must be transient, so that
#' \eqn{I_d-P_{\mathcal{D},\mathcal{D}}} is invertible. The asymptotic
#' variance additionally requires an irreducible chain and treats the initial
#' distribution as fixed.
#'
#' The initial probabilities on `down` are not renormalized. Paths starting
#' outside `down` have recovery time zero; to compute the mean conditional on
#' starting down, supply an initial distribution supported on `down`.
#'
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
#' @inherit var_mttr references
#'
#' @seealso [mttf()], [var_mttr()]
#'
#' @export
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
#' dimnames(P) <- list(c("working", "failed"), c("working", "failed"))
#' x <- MarkovChain(c(working = 0, failed = 1), P)
#' mttr(x, down = "failed")
#' var_mttr(x, down = "failed")
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
#' @param down Nonempty integer-valued numeric vector of state positions, or
#'   character vector of state names, giving the non-operational states.
#'   Entries must be distinct and select a proper subset of the states.
#'   Input order is ignored; states are used in their order in `x`.
#'
#' @details
#' These functions use [mttf()] and [var_mttf()] with `down` as the selected
#' state set. The non-operational subchain must be transient, so that
#' \eqn{I_d-P_{\mathcal{D},\mathcal{D}}} is invertible. The asymptotic
#' variance additionally requires an irreducible chain and treats the initial
#' distribution as fixed.
#'
#' The initial probabilities on `down` are not renormalized. Paths starting
#' outside `down` have recovery time zero; to compute the mean conditional on
#' starting down, supply an initial distribution supported on `down`.
#'
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
#' A numeric scalar containing the asymptotic variance \eqn{\sigma_\psi^2}
#' of the MTTR estimator. This is the variance in the \eqn{\sqrt{n}} limit;
#' the approximate variance of the estimator itself is \eqn{\sigma_\psi^2/n}.
#'
#' @references
#' Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026).
#' \emph{Matrix asymptotic calculus for plug-in maximum likelihood estimators in finite Markov chains}.
#' arXiv preprint arXiv:2607.17161. \url{https://arxiv.org/abs/2607.17161}.
#' See Proposition 11 for the MTTR differential and asymptotic variance.
#'
#' Sadek, A., and Limnios, N. (2002).
#' \emph{Asymptotic properties for maximum likelihood estimators for reliability and failure rates of Markov chains}.
#' Communications in Statistics - Theory and Methods, 31(10), 1837--1861.
#' \doi{10.1081/STA-120014916}.
#'
#' @seealso [mttr()], [var_mttf()]
#'
#' @export
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
#' dimnames(P) <- list(c("working", "failed"), c("working", "failed"))
#' x <- MarkovChain(c(working = 0, failed = 1), P)
#' mttr(x, down = "failed")
#' var_mttr(x, down = "failed")
var_mttr <- function(x, down) {
  var_mttf(x, down)
}
