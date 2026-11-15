# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Discrete-Time Intensity of Hitting Down                                    ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Discrete-Time Intensity of Hitting Down
#'
#' @description
#' Compute the probability of an operational-to-down transition at a given
#' positive time in a finite-state Markov chain.
#' @inheritParams ftp
#'
#' @details
#' Equation (48) of Gavrilopoulos, Trevezas, and Votsi (2026) gives
#' \deqn{r(t) = \Pr(X_{t-1}\in U, X_t\in D)
#'   = a P^{t-1} J_U P e_D, \qquad t\geq1,}
#' where \eqn{U} is the operational set, \eqn{D} its complement,
#' \eqn{J_U} is the diagonal operational-state indicator matrix, and
#' \eqn{e_D} is the down-state indicator vector.
#' The full transition matrix evolves the initial distribution to time
#' \eqn{t-1}; paths may have failed and recovered before that time.
#' This differs from the first-entry probability returned by [ftp()].
#' Paths starting down can contribute after recovery. The initial
#' distribution is fixed and is not renormalized.
#'
#' The vector \eqn{J_U P e_D} is computed by summing down-transition
#' probabilities directly and setting its down-state entries to zero.
#' Powers are computed once by [Ppowers()]. Empty and full operational
#' sets give zero intensity.
#'
#' @return A numeric vector of the same length and order as `time` containing
#'   the discrete-time intensities.
#' @inherit vcov_rel references
#' @seealso [ftp()], [avail()], [vcov_dti()]
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
#' x <- MarkovChain(c(0.4, 0.6), P)
#' dti(x, time = 1:3, up = 1L)
#' vcov_dti(x, time = 3L, up = 1L)
#' @export
dti <- function(x, time, up) {
  x <- validate_MarkovChain(x)
  validate_ftp_time(time)
  up <- validate_avail_up(up, get_state_names(x))
  down <- setdiff(get_state_names(x), up)
  if (length(up) == 0L || length(down) == 0L) return(numeric(length(time)))

  powers <- Ppowers(x$P, max(1L, time - 1L))
  exit <- setNames(numeric(nstates(x)), get_state_names(x))
  exit[up] <- rowSums(x$P[up, down, drop = FALSE])
  vapply(time, function(t) {
    as.numeric(x$a %*% powers[[t]] %*% exit)
  }, numeric(1), USE.NAMES = FALSE)
}

#' Asymptotic Variance of the Discrete-Time Intensity
#'
#' @description
#' Compute the pointwise asymptotic variance of the plug-in discrete-time
#' intensity estimator for an irreducible finite-state Markov chain,
#' treating the initial distribution as fixed.
#' @inheritParams vcov_ftp
#'
#' @details
#' Equation (59) in Proposition 12 of Gavrilopoulos, Trevezas, and Votsi
#' (2026) applies the differential
#' \deqn{L_t(H) = a\left\{
#'   \sum_{i=0}^{t-2} P^i H P^{t-2-i} J_U P
#'   + P^{t-1} J_U H\right\}e_D.}
#' The sum is empty at \eqn{t=1}. Both the distribution at time \eqn{t-1}
#' and the final transition contribute; their cross-covariance is included.
#' At later times, transitions from down states can affect the variance
#' through the first term.
#'
#' With \eqn{G} the matrix of entry-wise coefficients of \eqn{L_t}, let
#' \eqn{g=\operatorname{vec}(G^{\mathsf T})^{\mathsf T}} denote the
#' row-wise derivative. The variance \eqn{g\Sigma_{\mathbf p}g^{\mathsf T}}
#' is computed using the equivalent full-chain row-covariance contraction
#' \deqn{\sum_u \frac{1}{\pi_u}\sum_{j<l}
#'   p_{uj}p_{ul}(G_{uj}-G_{ul})^2,}
#' where \eqn{\pi} is the stationary distribution of the full chain.
#' This avoids constructing a dense transition covariance or subtracting
#' nearly equal quadratic terms. The derivative is assembled using outer
#' products, without a full matrix-power Jacobian.
#'
#' @return A numeric scalar containing the variance in the \eqn{\sqrt{n}}
#'   limit of the intensity estimator. Divide by the number of observed
#'   transitions \eqn{n} to approximate the variance of the estimator itself.
#'   This is a pointwise variance, not a joint covariance across times.
#' @inherit vcov_rel references
#' @seealso [dti()], [vcov_ftp()], [vcov_p()]
#' @inherit dti examples
#' @export
vcov_dti <- function(x, time, up) {
  x <- validate_MarkovChain(x)
  validate_ftp_time(time, scalar = TRUE)
  up <- validate_avail_up(up, get_state_names(x))
  down <- setdiff(get_state_names(x), up)
  stationary <- stationary_for_covariance(x$P)
  if (length(up) == 0L || length(down) == 0L) return(0)

  s <- nstates(x)
  state_names <- get_state_names(x)
  powers <- Ppowers(x$P, max(1L, time - 1L))
  exit <- setNames(numeric(s), state_names)
  exit[up] <- rowSums(x$P[up, down, drop = FALSE])
  G <- matrix(0, s, s, dimnames = list(state_names, state_names))
  # The last-transition differential a P^(t-1) J_U H e_D.
  previous <- as.numeric(x$a %*% powers[[time]])
  G[up, down] <- previous[match(up, state_names)]
  # The distribution differential a sum P^i H P^(t-2-i) J_U P e_D.
  for (i in seq_len(time - 1L)) {
    G <- G + outer(as.numeric(x$a %*% powers[[i]]),
                   as.numeric(powers[[time - i]] %*% exit))
  }

  variance <- 0
  for (j in seq_len(s - 1L)) {
    for (l in seq.int(j + 1L, s)) {
      variance <- variance + sum(x$P[, j] * x$P[, l] *
        (G[, j] - G[, l])^2 / stationary)
    }
  }
  as.numeric(variance)
}

#' Asymptotic Confidence Intervals for Discrete-Time Intensity
#'
#' Compute two-sided, pointwise normal-approximation confidence intervals using [dti()]
#' and [vcov_dti()].
#' @inheritParams dti
#' @inheritParams ci_rel
#' @inheritSection ci_rel Normal approximation
#' @details
#' The asymptotic normality result is Proposition 12 of Gavrilopoulos, Trevezas, and Votsi (2026).
#' The domain restrictions of [vcov_dti()] apply.
#' @return A data frame with numeric columns `estimate`, `std.error`,
#'   `lower`, and `upper`, giving the point estimate, estimated standard
#'   error, and confidence limits.
#'   There is one row per requested time, preserving order and duplicates.
#' @inherit ci_rel references
#' @seealso [dti()], [vcov_dti()], [ci_rel()]
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE)
#' x <- MarkovChain(c(0.4, 0.6), P)
#' ci_dti(x, time = 1:3, up = 1L, n = 1000)
#' @export
ci_dti <- function(x, time, up, n, level = 0.95) {
  validate_ci_parameters(n, level)
  estimate <- dti(x, time, up)
  times <- unique(time)
  variance <- vapply(times, function(t) vcov_dti(x, t, up),
                     numeric(1), USE.NAMES = FALSE)
  norm_ci(estimate, variance[match(time, times)], n, level)
}
