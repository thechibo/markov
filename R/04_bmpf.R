# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BMP Failure Rate                                                          ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' BMP Failure Rate of a Markov Chain
#'
#' @description
#' Compute the discrete-time BMP failure rate: the probability of first
#' failure at a given time conditional on survival to the preceding time.
#' @inheritParams rel
#'
#' @details
#' Equation (2.7) of Sadek and Limnios (2002) defines
#' \deqn{\lambda(t)=1-R(t)/R(t-1), \qquad t\geq1,}
#' when \eqn{R(t-1)>0}, and sets \eqn{\lambda(t)=0} otherwise.
#' At time zero, \eqn{\lambda(0)=1-R(0)} is the initial probability of
#' being down. Here \eqn{R(t)=a_U P_{UU}^t\mathbf{1}}.
#'
#' The implementation propagates the conditional distribution among surviving
#' operational states, normalizing at each step. Its weighted exit probability
#' equals the reliability ratio above, without subtracting nearly equal
#' probabilities or forming powers that may underflow. The initial distribution
#' is treated as fixed. Empty operational sets give one at time zero and zero
#' thereafter; full operational sets give zero at every time.
#'
#' @return A numeric vector of the same length and order as `time` containing
#'   the BMP failure rates, including the paper's zero-survival convention.
#' @inherit rel references
#' @seealso [rel()], [ftp()], [vcov_bmp()]
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
#' x <- MarkovChain(c(0.4, 0.6), P)
#' bmpf(x, time = 0:3, up = 1L)
#' vcov_bmp(x, time = 3L, up = 1L)
#' @export
bmpf <- function(x, time, up) {
  x <- validate_MarkovChain(x)
  validate_rel_time(time)
  up <- validate_avail_up(up, get_state_names(x))
  bmp_survival(x, time, up)$rate
}

#' Asymptotic Variance of the BMP Failure Rate
#'
#' @description
#' Compute the pointwise asymptotic variance of the plug-in BMP failure-rate
#' estimator for an irreducible finite-state Markov chain with fixed initial
#' distribution.
#' @inheritParams vcov_rel
#'
#' @details
#' Theorem 7 of Sadek and Limnios (2002, p. 1851) gives
#' \eqn{\sigma_\lambda^2(t)=g_t\Sigma_U g_t^{\mathsf T}}, where
#' \eqn{\Sigma_U} is the full-chain transition covariance restricted to
#' operational-to-operational coordinates in row-wise order. Equation (3.11)
#' gives the derivative
#' \deqn{g_t = -\frac{R(t-1)\,dR(t)-R(t)\,dR(t-1)}{R(t-1)^2}.}
#' Both reliability derivatives use the same transition matrix; their
#' cross-covariance is included. The implementation differentiates the
#' normalized survival recursion, an algebraically equivalent quotient rule
#' that avoids underflow of the reliability probabilities.
#'
#' With \eqn{G} the row-wise reshape of \eqn{g_t}, \eqn{Q=P_{UU}},
#' \eqn{b_i=\sum_{j\notin U}p_{ij}}, and full-chain stationary distribution
#' \eqn{\pi}, the covariance is contracted as
#' \deqn{\sum_{i\in U}\frac{1}{\pi_i}\left\{
#'   \sum_{j<l;\ j,l\in U}Q_{ij}Q_{il}(G_{ij}-G_{il})^2
#'   + b_i\sum_{j\in U}Q_{ij}G_{ij}^2\right\}.}
#' This avoids a dense transition covariance and subtractive cancellation.
#'
#' At positive times Theorem 7 requires \eqn{R(t-1)>0}; an error is raised
#' when this condition fails, even though [bmpf()] defines the rate as zero
#' there. At time zero the variance is zero because the initial distribution
#' is fixed. Full operational sets also have zero variance.
#'
#' @return A numeric scalar containing the variance in the \eqn{\sqrt{n}}
#'   limit of the BMP estimator. Divide by the number of observed transitions
#'   \eqn{n} to approximate the estimator variance. This is a pointwise
#'   variance, not a joint covariance across times.
#' @inherit rel references
#' @seealso [bmpf()], [vcov_rel()], [vcov_p()]
#' @inherit bmpf examples
#' @export
vcov_bmp <- function(x, time, up) {
  x <- validate_MarkovChain(x)
  validate_rel_time(time, scalar = TRUE)
  up <- validate_avail_up(up, get_state_names(x))
  stationary <- stationary_for_covariance(x$P)
  if (time == 0L || length(up) == nstates(x)) return(0)

  result <- bmp_survival(x, time, up, derivative = TRUE)
  r <- length(up)
  G <- matrix(result$gradient, r, r, byrow = TRUE)
  Q <- x$P[up, up, drop = FALSE]
  exit <- rowSums(x$P[up, setdiff(get_state_names(x), up), drop = FALSE])
  variance <- sum(exit * rowSums(Q * G^2) / stationary[up])
  for (j in seq_len(r - 1L)) {
    for (l in seq.int(j + 1L, r)) {
      variance <- variance + sum(Q[, j] * Q[, l] *
        (G[, j] - G[, l])^2 / stationary[up])
    }
  }
  as.numeric(variance)
}

# Conditional survival and its derivative with respect to row-wise Q entries.
bmp_survival <- function(x, time, up, derivative = FALSE) {
  r <- length(up)
  mass <- sum(x$a[up])
  rate <- numeric(length(time))
  survival <- rep(NA_real_, length(time))
  down <- setdiff(get_state_names(x), up)
  rate[time == 0L] <- sum(x$a[down])
  survival[time == 0L] <- mass
  if (length(down) == 0L) survival[] <- 1
  if (max(time) == 0L || length(down) == 0L) {
    return(list(rate = rate, survival = survival))
  }
  Q <- x$P[up, up, drop = FALSE]
  exit <- rowSums(x$P[up, down, drop = FALSE])
  v <- if (mass > 0) as.numeric(x$a[up]) / mass else numeric(r)
  if (derivative) J <- matrix(0, r, r^2)
  targets <- sort(unique(time[time > 0L]))
  step <- 1L
  for (target in targets) {
    while (step < target && mass > 0) {
      z <- as.numeric(v %*% Q)
      mass <- sum(z)
      if (mass > 0) {
        if (derivative) {
          B <- t(Q) %*% J
          for (u in seq_len(r)) {
            columns <- (u - 1L) * r + seq_len(r)
            B[, columns] <- B[, columns, drop = FALSE] + diag(v[u], r)
          }
          J <- (B - outer(z / mass, colSums(B))) / mass
        }
        v <- z / mass
      }
      step <- step + 1L
    }
    if (mass == 0) {
      if (derivative) {
        stop("BMP variance requires positive survival probability at `time - 1`.",
             call. = FALSE)
      }
      break
    }
    rate[time == target] <- sum(v * exit)
    survival[time == target] <- sum(v %*% Q)
  }
  if (!derivative) return(list(rate = rate, survival = survival))
  # sum(dv) = 0, so -rowSums(Q) %*% J equals exit %*% J.
  gradient <- as.numeric(exit %*% J) - rep(v, each = r)
  list(rate = rate, survival = survival, gradient = gradient)
}

#' Asymptotic Confidence Intervals for BMP Failure Rate
#'
#' Compute two-sided, pointwise normal-approximation confidence intervals using [bmpf()]
#' and [vcov_bmp()].
#' @inheritParams bmpf
#' @inheritParams ci_rel
#' @inheritSection ci_rel Normal approximation
#' @details
#' The asymptotic normality result is Theorem 7 of Sadek and Limnios (2002).
#' The domain restrictions of [vcov_bmp()] apply.
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
#' @seealso [bmpf()], [vcov_bmp()], [ci_rel()]
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE)
#' x <- MarkovChain(c(0.4, 0.6), P)
#' ci_bmpf(x, time = 1:3, up = 1L, n = 1000)
#' @export
ci_bmpf <- function(x, time, up, n, level = 0.95) {
  validate_ci_parameters(n, level)
  estimate <- bmpf(x, time, up)
  times <- unique(time)
  variance <- vapply(times, function(t) vcov_bmp(x, t, up),
                     numeric(1), USE.NAMES = FALSE)
  norm_ci(estimate, variance[match(time, times)], n, level)
}
