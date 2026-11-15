# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stationary Availability                                                   ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Stationary Availability of a Markov Chain
#'
#' @description
#' Compute the stationary probability of being operational for an
#' irreducible finite-state Markov chain.
#' @inheritParams rel
#'
#' @details
#' Equation (49) of Gavrilopoulos, Trevezas, and Votsi (2026) gives
#' \deqn{A_\infty = \pi e_U = \sum_{i\in U}\pi_i,}
#' where \eqn{\pi} is the stationary distribution of the full transition
#' matrix and \eqn{e_U} is the operational-state indicator vector.
#' This function uses the same stationary solver as [stat()]. Both this
#' function and [vcov_savail()] require an irreducible chain, so the
#' stationary distribution is unique. The result does not depend on the
#' initial distribution, although the complete `MarkovChain` is validated.
#'
#' Periodic chains are supported. For them, stationary availability need
#' not equal the limit of [avail()] from an arbitrary initial distribution.
#' Empty and full operational sets give zero and one, respectively.
#'
#' @return A numeric scalar containing the stationary availability.
#' @inherit vcov_rel references
#' @seealso [stat()], [avail()], [vcov_savail()]
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
#' x <- MarkovChain(c(0.4, 0.6), P)
#' savail(x, up = 1L)
#' vcov_savail(x, up = 1L)
#' @export
savail <- function(x, up) {
  x <- validate_MarkovChain(x)
  up <- validate_avail_up(up, get_state_names(x))
  stationary <- stationary_for_covariance(x$P)
  if (length(up) == nstates(x)) return(1)
  as.numeric(sum(stationary[up]))
}

#' Asymptotic Variance of Stationary Availability
#'
#' @description
#' Compute the asymptotic variance of the plug-in stationary availability
#' estimator for an irreducible finite-state Markov chain.
#' @inheritParams savail
#'
#' @details
#' Equation (60) in Proposition 12 of Gavrilopoulos, Trevezas, and Votsi
#' (2026) gives the limiting variable \eqn{W_\pi e_U}. Consequently,
#' \deqn{\sigma^2_{A_\infty} = e_U^{\mathsf T}\Sigma_\pi e_U,}
#' where \eqn{\Sigma_\pi} is the stationary-distribution asymptotic
#' covariance from Proposition 7, computed by [vcov_stat()]. The function
#' sums all entries of its operational-state submatrix, including the
#' off-diagonal covariances. Since the stationary probabilities sum to one,
#' the complementary state set has the same variance; the smaller set is
#' used to reduce summation and cancellation. Empty and full sets return
#' exactly zero after checking irreducibility.
#'
#' The result does not depend on the initial distribution. This is the
#' asymptotic variance of the transition-MLE plug-in estimator of stationary
#' availability. The numerical conditioning limitations of [vcov_stat()]
#' also apply here.
#'
#' @return A numeric scalar containing the variance in the \eqn{\sqrt{n}}
#'   limit of the stationary availability estimator. Divide by the number
#'   of observed transitions \eqn{n} to approximate the variance of the
#'   estimator itself.
#' @inherit vcov_rel references
#' @seealso [savail()], [vcov_stat()], [vcov_avail()]
#' @inherit savail examples
#' @export
vcov_savail <- function(x, up) {
  x <- validate_MarkovChain(x)
  state_names <- get_state_names(x)
  up <- validate_avail_up(up, state_names)
  if (length(up) == 0L || length(up) == nstates(x)) {
    stationary_for_covariance(x$P)
    return(0)
  }

  if (length(up) > nstates(x) / 2) up <- setdiff(state_names, up)
  positions <- match(up, state_names)
  Sigma <- vcov_stat(x)
  as.numeric(sum(Sigma[positions, positions, drop = FALSE]))
}

#' Asymptotic Confidence Intervals for Stationary Availability
#'
#' Compute two-sided, pointwise normal-approximation confidence intervals using [savail()]
#' and [vcov_savail()].
#' @inheritParams savail
#' @inheritParams ci_rel
#' @inheritSection ci_rel Normal approximation
#' @details
#' The asymptotic normality result is Proposition 12 of Gavrilopoulos, Trevezas, and Votsi (2026).
#' The domain restrictions of [vcov_savail()] apply.
#' @return A data frame with numeric columns `estimate`, `std.error`,
#'   `lower`, and `upper`, giving the point estimate, estimated standard
#'   error, and confidence limits.
#'   There is one row.
#' @inherit ci_rel references
#' @seealso [savail()], [vcov_savail()], [ci_rel()]
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE)
#' x <- MarkovChain(c(0.4, 0.6), P)
#' ci_savail(x, up = 1L, n = 1000)
#' @export
ci_savail <- function(x, up, n, level = 0.95) {
  validate_ci_parameters(n, level)
  estimate <- savail(x, up)
  norm_ci(estimate, vcov_savail(x, up), n, level)
}
