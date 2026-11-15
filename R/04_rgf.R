# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RG Failure Rate                                                           ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' RG Failure Rate of a Markov Chain
#'
#' @description
#' Compute the discrete-time RG failure rate as the logarithmic decrease
#' in reliability over one time step.
#' @inheritParams rel
#'
#' @details
#' Equation (2.9) of Sadek and Limnios (2002, p. 1842) defines
#' \deqn{r(t)=\log\{R(t-1)/R(t)\}, \quad t\geq1, \qquad
#'   r(0)=-\log R(0).}
#' When the reliability ratio is defined, this equals
#' \eqn{-\log(1-\lambda(t))}, with BMP rate \eqn{\lambda(t)}.
#' The implementation reuses the normalized survival recursion underlying
#' [bmpf()], avoiding underflow of the raw reliability powers. It uses
#' `log1p()` for small failure probabilities and the directly computed
#' conditional survival probability when the BMP rate is close to one.
#' The initial distribution is treated as fixed.
#'
#' If \eqn{R(t-1)>0} but \eqn{R(t)=0}, the rate is `Inf`. If both are
#' zero, the defining ratio is undefined and `NaN` is returned, rather
#' than applying the zero-rate convention of [bmpf()]. At time zero,
#' \eqn{R(0)=0} gives `Inf`. Full operational sets have rate zero.
#'
#' @return A numeric vector of the same length and order as `time`,
#'   containing RG rates, with `Inf` and `NaN` as described above.
#' @inherit rel references
#' @seealso [bmpf()], [rel()], [vcov_rg()]
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
#' x <- MarkovChain(c(0.4, 0.6), P)
#' rgf(x, time = 0:3, up = 1L)
#' vcov_rg(x, time = 3L, up = 1L)
#' @export
rgf <- function(x, time, up) {
  x <- validate_MarkovChain(x)
  validate_rel_time(time)
  up <- validate_avail_up(up, get_state_names(x))
  profile <- bmp_survival(x, time, up)
  value <- rep(NaN, length(time))
  defined <- !is.na(profile$survival)
  small <- defined & profile$rate <= 0.5
  value[small] <- -log1p(-profile$rate[small])
  value[defined & !small] <- -log(profile$survival[defined & !small])
  value
}

#' Asymptotic Variance of the RG Failure Rate
#'
#' @description
#' Compute the pointwise asymptotic variance of the plug-in RG failure-rate
#' estimator for an irreducible finite-state Markov chain with fixed initial
#' distribution.
#' @inheritParams vcov_rel
#'
#' @details
#' Corollary 2 of Sadek and Limnios (2002, pp. 1851--1852) applies the delta
#' method to \eqn{r(t)=-\log(1-\lambda(t))}, giving
#' \deqn{\sigma_r^2(t)=\frac{\sigma_\lambda^2(t)}{(1-\lambda(t))^2}.}
#' This function reuses [vcov_bmp()] for \eqn{\sigma_\lambda^2(t)} and the
#' normalized survival recursion for \eqn{1-\lambda(t)=R(t)/R(t-1)}.
#' The scaling is performed by two successive divisions, avoiding premature
#' underflow from squaring a small survival probability.
#'
#' At positive times, both preceding and current survival must be positive.
#' An error is raised for infinite or undefined RG rates, where Corollary 2
#' does not provide a finite-rate delta-method limit. At time zero the
#' variance is zero when \eqn{R(0)>0}, because the initial distribution is
#' fixed. Full operational sets also have variance zero. Very small survival
#' ratios can produce an asymptotic variance exceeding the numeric range.
#'
#' @return A numeric scalar containing the variance in the \eqn{\sqrt{n}}
#'   limit of the RG estimator. Divide by the number of observed transitions
#'   \eqn{n} to approximate the estimator variance. This is a pointwise
#'   variance, not a joint covariance across times.
#' @inherit rel references
#' @seealso [rgf()], [vcov_bmp()]
#' @inherit rgf examples
#' @export
vcov_rg <- function(x, time, up) {
  x <- validate_MarkovChain(x)
  validate_rel_time(time, scalar = TRUE)
  up <- validate_avail_up(up, get_state_names(x))
  survival <- bmp_survival(x, time, up)$survival
  if (is.na(survival) || survival <= 0) {
    stop("RG variance requires positive survival and a finite RG rate.",
         call. = FALSE)
  }
  variance <- vcov_bmp(x, time = time, up = up)
  as.numeric(variance / survival / survival)
}

#' Asymptotic Confidence Intervals for RG Failure Rate
#'
#' Compute two-sided, pointwise normal-approximation confidence intervals using [rgf()]
#' and [vcov_rg()].
#' @inheritParams rgf
#' @inheritParams ci_rel
#' @inheritSection ci_rel Normal approximation
#' @details
#' The asymptotic normality result is Corollary 2 of Sadek and Limnios (2002).
#' The domain restrictions of [vcov_rg()] apply.
#' @return A data frame with numeric columns `estimate`, `std.error`,
#'   `lower`, and `upper`, giving the point estimate, estimated standard
#'   error, and confidence limits.
#'   There is one row per requested time, preserving order and duplicates.
#' @references
#' Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026).
#' \emph{Matrix asymptotic calculus for plug-in maximum likelihood estimators in finite Markov chains}.
#' arXiv preprint arXiv:2607.17161v2. \url{https://arxiv.org/abs/2607.17161v2}.
#'
#' Sadek, A., and Limnios, N. (2002).
#' \emph{Asymptotic properties for maximum likelihood estimators for reliability and failure rates of Markov chains}.
#' Communications in Statistics - Theory and Methods, 31(10), 1837--1861.
#' \doi{10.1081/STA-120014916}.
#' @seealso [rgf()], [vcov_rg()], [ci_rel()]
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE)
#' x <- MarkovChain(c(0.4, 0.6), P)
#' ci_rgf(x, time = 1:3, up = 1L, n = 1000)
#' @export
ci_rgf <- function(x, time, up, n, level = 0.95) {
  validate_ci_parameters(n, level)
  estimate <- rgf(x, time, up)
  times <- unique(time)
  variance <- vapply(times, function(t) vcov_rg(x, t, up),
                     numeric(1), USE.NAMES = FALSE)
  norm_ci(estimate, variance[match(time, times)], n, level)
}
