# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Failure-Time Probability                                                  ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Failure-Time Probability of a Markov Chain
#'
#' @description
#' Compute the probability that a finite-state Markov chain first enters a
#' non-operational state at the requested positive time.
#' @inheritParams rel
#' @param time A nonempty positive integer vector of discrete time points,
#'   each less than `.Machine$integer.max`. Missing and infinite values are
#'   not allowed.
#'
#' @details
#' Let \eqn{U} be the operational states, \eqn{D} their complement, and
#' \eqn{T_D = \inf\{t \geq 0 : X_t \in D\}}. Equation (44) of
#' Gavrilopoulos, Trevezas, and Votsi (2026) gives
#' \deqn{f(t) = \Pr(T_D=t) = a_U Q^{t-1} b, \qquad
#'   Q=P_{UU}, \quad b=P_{UD}\mathbf{1}_{|D|}, \quad t\geq1.}
#' Neither the initial distribution nor the restricted block is renormalized.
#' Paths starting down have failure time zero and do not contribute at the
#' positive times supported here. Their probability is \eqn{\sum_{i\in D}a_i}.
#' The exit vector is summed directly from \eqn{P_{UD}}, avoiding subtraction
#' of nearly equal reliability probabilities. Powers are computed by
#' [Ppowers()]. Empty and full operational sets give zero at positive times.
#'
#' @return A numeric vector of the same length and order as `time` containing
#'   the failure-time probabilities.
#' @inherit vcov_rel references
#' @seealso [rel()], [vcov_ftp()]
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
#' x <- MarkovChain(c(0.4, 0.6), P)
#' ftp(x, time = 1:3, up = 1L)
#' vcov_ftp(x, time = 3L, up = 1L)
#' @export
ftp <- function(x, time, up) {
  x <- validate_MarkovChain(x)
  validate_ftp_time(time)
  up <- validate_avail_up(up, get_state_names(x))
  down <- setdiff(get_state_names(x), up)
  if (length(up) == 0L || length(down) == 0L) return(numeric(length(time)))

  powers <- Ppowers(x$P[up, up, drop = FALSE], max(1L, time - 1L))
  exit <- rowSums(x$P[up, down, drop = FALSE])
  vapply(time, function(t) {
    as.numeric(x$a[up] %*% powers[[t]] %*% exit)
  }, numeric(1), USE.NAMES = FALSE)
}

#' Asymptotic Variance of the Failure-Time Probability
#'
#' @description
#' Compute the pointwise asymptotic variance of the plug-in failure-time
#' probability estimator for an irreducible finite-state Markov chain,
#' treating the initial distribution as fixed.
#' @inheritParams ftp
#' @param time A single positive integer less than `.Machine$integer.max`.
#'   Missing and infinite values are not allowed.
#'
#' @details
#' Equation (53) in Proposition 10 of Gavrilopoulos, Trevezas, and Votsi
#' (2026) gives the limiting variable by applying the differential
#' \deqn{L_t(H) = a_U\left\{
#'   \sum_{i=0}^{t-2} Q^i H_{UU} Q^{t-2-i} b
#'   + Q^{t-1} H_{UD}\mathbf{1}_{|D|}\right\},}
#' where \eqn{Q=P_{UU}}, \eqn{b=P_{UD}\mathbf{1}_{|D|}}, and the sum is
#' empty at \eqn{t=1}. Both terms use the same transition estimation error;
#' their cross-covariance is included.
#'
#' If \eqn{G} contains the entry-wise coefficients of \eqn{L_t}, the
#' variance is \eqn{g\Sigma_{\mathbf p}g^{\mathsf T}}, with
#' \eqn{g=\operatorname{vec}(G^{\mathsf T})^{\mathsf T}}. The implementation
#' contracts the full-chain multinomial row covariances using the equivalent
#' expression
#' \deqn{\sum_{u\in U}\frac{1}{\pi_u}
#'   \sum_{j<l}p_{uj}p_{ul}(G_{uj}-G_{ul})^2.}
#' This includes cross-covariances between the operational and exit blocks,
#' avoids constructing a dense transition covariance, and avoids subtracting
#' nearly equal quadratic terms. The stationary probabilities \eqn{\pi_u}
#' refer to the full chain. No new matrix dependencies are required.
#'
#' @return A numeric scalar containing the variance in the \eqn{\sqrt{n}}
#'   limit of the failure-time probability estimator. Divide by the number
#'   of observed transitions \eqn{n} to approximate the variance of the
#'   estimator itself. This is a pointwise variance, not a joint covariance
#'   across times.
#' @inherit vcov_rel references
#' @seealso [ftp()], [vcov_rel()], [vcov_p()]
#' @inherit ftp examples
#' @export
vcov_ftp <- function(x, time, up) {
  x <- validate_MarkovChain(x)
  validate_ftp_time(time, scalar = TRUE)
  up <- validate_avail_up(up, get_state_names(x))
  down <- setdiff(get_state_names(x), up)
  stationary <- stationary_for_covariance(x$P)
  if (length(up) == 0L || length(down) == 0L) return(0)

  powers <- Ppowers(x$P[up, up, drop = FALSE], max(1L, time - 1L))
  exit <- rowSums(x$P[up, down, drop = FALSE])
  G <- matrix(0, nrow = length(up), ncol = nstates(x),
              dimnames = list(up, get_state_names(x)))
  # The exit-block differential a_U Q^(t-1) H_UD 1.
  G[, down] <- as.numeric(x$a[up] %*% powers[[time]])
  # The power differential a_U sum Q^i H_UU Q^(t-2-i) b.
  for (i in seq_len(time - 1L)) {
    G[, up] <- G[, up, drop = FALSE] +
      outer(as.numeric(x$a[up] %*% powers[[i]]),
            as.numeric(powers[[time - i]] %*% exit))
  }

  variance <- 0
  P <- x$P[up, , drop = FALSE]
  for (j in seq_len(nstates(x) - 1L)) {
    for (l in seq.int(j + 1L, nstates(x))) {
      variance <- variance + sum(P[, j] * P[, l] *
        (G[, j] - G[, l])^2 / stationary[up])
    }
  }
  as.numeric(variance)
}

validate_ftp_time <- function(time, scalar = FALSE) {
  validate_rel_time(time, scalar = scalar)
  if (any(time == 0L)) {
    stop("`time` must contain positive integers.", call. = FALSE)
  }
  invisible(NULL)
}

#' Asymptotic Confidence Intervals for Failure-Time Probability
#'
#' Compute two-sided, pointwise normal-approximation confidence intervals using [ftp()]
#' and [vcov_ftp()].
#' @inheritParams ftp
#' @inheritParams ci_rel
#' @inheritSection ci_rel Normal approximation
#' @details
#' The asymptotic normality result is Proposition 10 of Gavrilopoulos, Trevezas, and Votsi (2026).
#' The domain restrictions of [vcov_ftp()] apply.
#' @return A data frame with numeric columns `estimate`, `std.error`,
#'   `lower`, and `upper`, giving the point estimate, estimated standard
#'   error, and confidence limits.
#'   There is one row per requested time, preserving order and duplicates.
#' @inherit ci_rel references
#' @seealso [ftp()], [vcov_ftp()], [ci_rel()]
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE)
#' x <- MarkovChain(c(0.4, 0.6), P)
#' ci_ftp(x, time = 1:3, up = 1L, n = 1000)
#' @export
ci_ftp <- function(x, time, up, n, level = 0.95) {
  validate_ci_parameters(n, level)
  estimate <- ftp(x, time, up)
  times <- unique(time)
  variance <- vapply(times, function(t) vcov_ftp(x, t, up),
                     numeric(1), USE.NAMES = FALSE)
  norm_ci(estimate, variance[match(time, times)], n, level)
}
