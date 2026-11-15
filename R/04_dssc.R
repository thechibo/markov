# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Down-State Survival Curve                                                 ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Down-State Survival Curve of a Markov Chain
#'
#' @description
#' Compute the probability that a finite-state Markov chain remains in the
#' non-operational states from time zero through the requested time.
#' @inheritParams rel
#' @param down Integer-valued numeric vector of state positions or character
#'   vector of state names giving the non-operational states. Missing values
#'   are not allowed. Duplicates are ignored; empty and full state sets are
#'   valid. States are used in their order in `x`.
#'
#' @details
#' For the set of down states \eqn{D}, the down-state survival curve is
#' \deqn{Q(t) = \Pr\{X_0,\ldots,X_t \in D\}
#'   = a_D P_{DD}^t \mathbf{1}_d, \qquad d = |D|.}
#' This is the reliability formula in Equation (2.5) of Sadek and Limnios
#' (2002), with the selected states changed from operational to
#' non-operational. The function calls [rel()] with `up = down`.
#' Neither the restricted initial distribution nor the transition block is
#' renormalized. Thus \eqn{Q(0) = \sum_{i \in D} a_i}; to obtain survival
#' conditional on starting down, supply an initial distribution supported
#' on `down`. Paths that leave `down` and later return do not contribute.
#' The curve is not generally the complement of reliability on the up states.
#'
#' @return A numeric vector of the same length and order as `time` containing
#'   the down-state survival probabilities.
#' @inherit rel references
#' @seealso [rel()], [mttr()], [vcov_dssc()], [vcov_dssc_sl()]
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
#' x <- MarkovChain(c(0.4, 0.6), P)
#' dssc(x, time = 0:3, down = 2L)
#' vcov_dssc(x, time = 3L, down = 2L)
#' vcov_dssc_sl(x, time = 3L, down = 2L)
#' @export
dssc <- function(x, time, down) {
  rel(x, time = time, up = down)
}

#' Asymptotic Variance of Down-State Survival
#'
#' @description
#' Compute the pointwise asymptotic variance of the plug-in down-state
#' survival estimator for an irreducible finite-state Markov chain,
#' treating the initial distribution as fixed.
#' @inheritParams dssc
#' @param time A single nonnegative integer giving the discrete time point.
#'   Values must be less than `.Machine$integer.max`; missing and infinite
#'   values are not allowed.
#'
#' @details
#' Proposition 10, Equation (54), of Gavrilopoulos, Trevezas, and Votsi
#' (2026) gives the down-state survival limit. Its variance is Equation (52)
#' with the operational block replaced by the down-state block:
#' \deqn{\sigma_Q^2(t) = D_t \Sigma_D D_t^{\mathsf T}, \qquad
#'   D_t = \sum_{i=0}^{t-1} (a_D P_{DD}^i) \otimes
#'   (P_{DD}^{t-1-i}\mathbf{1}_d)^{\mathsf T}.}
#' The covariance \eqn{\Sigma_D} selects the down-to-down transition
#' coordinates in row-wise order from the full-chain transition covariance.
#' The function calls [vcov_rel()] with `up = down`, retaining the full-chain
#' stationary distribution and the unnormalized restricted initial
#' distribution. The variance at time zero is zero.
#'
#' @return A numeric scalar containing the variance in the \eqn{\sqrt{n}}
#'   limit of the down-state survival estimator. Divide by the number of
#'   observed transitions \eqn{n} to approximate the variance of the estimator
#'   itself. This is a pointwise variance, not a joint covariance across times.
#' @inherit vcov_rel references
#' @seealso [dssc()], [vcov_rel()], [vcov_dssc_sl()]
#' @inherit dssc examples
#' @export
vcov_dssc <- function(x, time, down) {
  vcov_rel(x, time = time, up = down)
}

#' Component-wise Reference for Down-State Survival Variance
#'
#' @description
#' Compute the pointwise asymptotic variance of down-state survival using
#' the component-wise covariance expansion of Sadek and Limnios (2002).
#' The full Markov chain must be irreducible.
#' @inheritParams vcov_dssc
#'
#' @details
#' This function calls [vcov_rel_sl()] with `up = down`, applying Theorem 6
#' (p. 1850) to the down-state block. It sums the weighted diagonal variances
#' and twice every unordered pairwise covariance of the entries of
#' \eqn{P_{DD}^t}. As in [vcov_rel_sl()], all unordered pairs are included
#' to implement the variance of the weighted sum; the printed simultaneous
#' coordinate-wise pair bounds would omit crossed pairs.
#' The initial distribution is fixed and its restriction to `down` is not
#' renormalized. Time zero has variance zero. This reference constructs a
#' dense component covariance and is intended for small down-state blocks.
#'
#' @inherit vcov_dssc return
#' @inherit rel references
#' @seealso [dssc()], [vcov_dssc()], [vcov_rel_sl()]
#' @inherit dssc examples
#' @export
vcov_dssc_sl <- function(x, time, down) {
  vcov_rel_sl(x, time = time, up = down)
}

#' Asymptotic Confidence Intervals for Down-State Survival
#'
#' Compute two-sided, pointwise normal-approximation confidence intervals using [dssc()]
#' and [vcov_dssc()].
#' @inheritParams dssc
#' @inheritParams ci_rel
#' @inheritSection ci_rel Normal approximation
#' @details
#' The asymptotic normality result is Proposition 10 of Gavrilopoulos, Trevezas, and Votsi (2026).
#' The domain restrictions of [vcov_dssc()] apply.
#' @return A data frame with numeric columns `estimate`, `std.error`,
#'   `lower`, and `upper`, giving the point estimate, estimated standard
#'   error, and confidence limits.
#'   There is one row per requested time, preserving order and duplicates.
#' @inherit ci_rel references
#' @seealso [dssc()], [vcov_dssc()], [ci_rel()]
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE)
#' x <- MarkovChain(c(0.4, 0.6), P)
#' ci_dssc(x, time = 1:3, down = 1L, n = 1000)
#' @export
ci_dssc <- function(x, time, down, n, level = 0.95) {
  validate_ci_parameters(n, level)
  estimate <- dssc(x, time, down)
  times <- unique(time)
  variance <- vapply(times, function(t) vcov_dssc(x, t, down),
                     numeric(1), USE.NAMES = FALSE)
  norm_ci(estimate, variance[match(time, times)], n, level)
}
