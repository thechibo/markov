# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Maintainability                                                           ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Maintainability of a Markov Chain
#'
#' @description
#' Compute the probability that a finite-state Markov chain has reached an
#' operational state by the requested time, including time zero.
#' @inheritParams dssc
#'
#' @details
#' Let \eqn{D} be the set of non-operational states and
#' \eqn{T_U = \inf\{t \geq 0 : X_t \notin D\}}. Maintainability is
#' \deqn{M(t) = \Pr(T_U \leq t) = 1 - Q(t),}
#' where \eqn{Q(t)} is the down-state survival curve. This function computes
#' `1 - dssc(x, time, down)`, as in Proposition 10 of Gavrilopoulos,
#' Trevezas, and Votsi (2026).
#'
#' The initial distribution is not renormalized. Consequently,
#' \eqn{M(0) = 1 - \sum_{i \in D} a_i}: paths starting operational have
#' recovery time zero. To obtain maintainability conditional on starting
#' down, supply an initial distribution supported on `down`. Paths that
#' recover and subsequently fail still contribute to maintainability.
#'
#' @return A numeric vector of the same length and order as `time` containing
#'   the maintainability probabilities.
#' @inherit vcov_rel references
#' @seealso [dssc()], [mttr()], [vcov_main()], [vcov_main_sl()]
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
#' x <- MarkovChain(c(0, 1), P)
#' main(x, time = 0:3, down = 2L)
#' vcov_main(x, time = 3L, down = 2L)
#' vcov_main_sl(x, time = 3L, down = 2L)
#' @export
main <- function(x, time, down) {
  1 - dssc(x, time = time, down = down)
}

#' Asymptotic Variance of Maintainability
#'
#' @description
#' Compute the pointwise asymptotic variance of the plug-in maintainability
#' estimator for an irreducible finite-state Markov chain, treating the
#' initial distribution as fixed.
#' @inheritParams vcov_dssc
#'
#' @details
#' Since \eqn{\widehat M(t) - M(t) = -(\widehat Q(t) - Q(t))}, the
#' asymptotic variance of maintainability equals that of down-state survival.
#' Proposition 10 of Gavrilopoulos, Trevezas, and Votsi (2026) states this
#' sign change after Equation (54). The variance is obtained from Equation
#' (52) with the down-state block substituted for the operational block.
#' This function calls [vcov_dssc()] without changing its result. The
#' full-chain stationary distribution is used and the variance at time zero
#' is zero.
#'
#' @return A numeric scalar containing the variance in the \eqn{\sqrt{n}}
#'   limit of the maintainability estimator. Divide by the number of observed
#'   transitions \eqn{n} to approximate the variance of the estimator itself.
#'   This is a pointwise variance, not a joint covariance across times.
#' @inherit vcov_rel references
#' @seealso [main()], [vcov_dssc()], [vcov_main_sl()]
#' @inherit main examples
#' @export
vcov_main <- function(x, time, down) {
  vcov_dssc(x, time = time, down = down)
}

#' Component-wise Reference for Maintainability Variance
#'
#' @description
#' Compute the pointwise asymptotic variance of maintainability using the
#' component-wise covariance expansion of Sadek and Limnios (2002).
#' The full Markov chain must be irreducible.
#' @inheritParams vcov_dssc
#'
#' @details
#' This function calls [vcov_dssc_sl()] because subtracting down-state
#' survival from one leaves its variance unchanged. The underlying reference
#' applies Theorem 6 (p. 1850) of Sadek and Limnios (2002) to the down-state
#' block, including every unordered pairwise covariance as documented in
#' [vcov_rel_sl()]. The initial distribution is fixed and is not renormalized
#' on `down`. Time zero has variance zero. The reference constructs a dense
#' component covariance and is intended for small down-state blocks.
#'
#' @inherit vcov_main return
#' @inherit rel references
#' @seealso [main()], [vcov_main()], [vcov_dssc_sl()]
#' @inherit main examples
#' @export
vcov_main_sl <- function(x, time, down) {
  vcov_dssc_sl(x, time = time, down = down)
}
