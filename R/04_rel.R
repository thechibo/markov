# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reliability                                                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Reliability of a Markov Chain
#'
#' @description
#' Compute the probability that a finite-state Markov chain remains in the
#' operational states from time zero through the requested time.
#'
#' @param x An object of class \code{MarkovChain}.
#' @param time A nonempty nonnegative integer vector of discrete time points.
#'   Values must be less than `.Machine$integer.max`; missing and infinite
#'   values are not allowed.
#' @param up Integer-valued numeric vector of state positions or character
#'   vector of state names giving the operational states. Missing values are
#'   not allowed. Duplicates are ignored; empty and full state sets are valid.
#'   States are used in their order in `x`.
#'
#' @details
#' Equation (2.5) of Sadek and Limnios (2002) gives
#' \deqn{R(t) = a_U Q^t \mathbf{1}_r, \qquad Q = P_{UU},}
#' where \eqn{a_U} is the initial distribution restricted to the operational
#' states and \eqn{r} is their number. Neither \eqn{a_U} nor \eqn{Q} is
#' renormalized. In particular, \eqn{R(0) = \sum_{i \in U} a_i}.
#' Paths that leave `up` and subsequently return do not contribute to
#' reliability. Powers of the restricted block are computed by [Ppowers()],
#' without an eigendecomposition.
#'
#' @return A numeric vector of the same length and order as `time` containing
#'   the reliability probabilities.
#' @references
#' Sadek, A., and Limnios, N. (2002).
#' \emph{Asymptotic properties for maximum likelihood estimators for reliability and failure rates of Markov chains}.
#' Communications in Statistics - Theory and Methods, 31(10), 1837--1861.
#' \doi{10.1081/STA-120014916}.
#' @seealso [avail()], [mttf()], [vcov_rel()], [vcov_rel_sl()]
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
#' x <- MarkovChain(c(0.4, 0.6), P)
#' rel(x, time = 0:3, up = 1L)
#' vcov_rel(x, time = 3L, up = 1L)
#' vcov_rel_sl(x, time = 3L, up = 1L)
#' @export
rel <- function(x, time, up) {
  x <- validate_MarkovChain(x)
  validate_rel_time(time)
  up <- validate_avail_up(up, get_state_names(x))

  if (length(up) == 0L) return(numeric(length(time)))
  if (length(up) == nstates(x)) return(rep(1, length(time)))

  Q <- x$P[up, up, drop = FALSE]
  powers <- Ppowers(Q, max(1L, time))
  ones <- rep(1, length(up))
  vapply(time, function(t) {
    as.numeric(x$a[up] %*% powers[[t + 1L]] %*% ones)
  }, numeric(1), USE.NAMES = FALSE)
}

#' Asymptotic Variance of Reliability
#'
#' @description
#' Compute the pointwise asymptotic variance of the plug-in reliability
#' estimator for an irreducible finite-state Markov chain, treating the
#' initial distribution as fixed.
#' @inheritParams rel
#' @param time A single nonnegative integer giving the discrete time point.
#'   Values must be less than `.Machine$integer.max`; missing and infinite
#'   values are not allowed.
#'
#' @details
#' Proposition 10, Equation (52), of Gavrilopoulos, Trevezas, and Votsi (2026)
#' gives the variance as \eqn{D_t \Sigma_U D_t^{\mathsf T}}, where
#' \deqn{D_t = \sum_{i=0}^{t-1}
#'   (a_U Q^i) \otimes (Q^{t-1-i}\mathbf{1}_r)^{\mathsf T}.}
#' Here \eqn{Q=P_{UU}} and \eqn{\Sigma_U} selects the operational-transition
#' coordinates from the full-chain transition covariance in row-wise order.
#' The Kronecker operator is contracted with the initial distribution and
#' vector of ones before forming the variance, avoiding a dense
#' \eqn{r^2 \times r^2} power derivative. The full-chain stationary
#' distribution is used; the operational block is not renormalized.
#' At time zero the derivative and variance are zero.
#'
#' @return A numeric scalar containing the variance in the \eqn{\sqrt{n}}
#'   limit of the reliability estimator. Divide by the number of observed
#'   transitions \eqn{n} to approximate the variance of the estimator itself.
#'   This is a pointwise variance, not a joint covariance across times.
#' @references
#' Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026).
#' \emph{Matrix asymptotic calculus for plug-in maximum likelihood estimators in finite Markov chains}.
#' arXiv preprint arXiv:2607.17161. \url{https://arxiv.org/abs/2607.17161}.
#' @seealso [rel()], [vcov_rel_sl()], [vcov_p()]
#' @inherit rel examples
#' @export
vcov_rel <- function(x, time, up) {
  x <- validate_MarkovChain(x)
  validate_rel_time(time, scalar = TRUE)
  up <- validate_avail_up(up, get_state_names(x))
  Sigma <- rel_transition_covariance(x, up)

  r <- length(up)
  if (time == 0L || r == 0L || r == nstates(x)) return(0)

  powers <- Ppowers(x$P[up, up, drop = FALSE], max(1L, time - 1L))
  D <- matrix(0, nrow = 1, ncol = r^2)
  ones <- rep(1, r)
  for (i in seq_len(time)) {
    D <- D + kp(x$a[up] %*% powers[[i]],
                t(powers[[time - i + 1L]] %*% ones))
  }
  as.numeric(D %*% Sigma %*% t(D))
}

#' Component-wise Reference for Reliability Variance
#'
#' @description
#' Compute the pointwise asymptotic variance using the component-wise
#' covariance expansion of Theorem 6 of Sadek and Limnios (2002, p. 1850).
#' The full Markov chain must be irreducible.
#' @inheritParams vcov_rel
#' @details
#' The covariance of the entries of \eqn{Q^t}, with \eqn{Q=P_{UU}}, is
#' computed using the component derivatives
#' \deqn{\frac{\partial (Q^t)_{ij}}{\partial Q_{uv}} =
#'   \sum_{h=0}^{t-1} (Q^h)_{iu}(Q^{t-1-h})_{vj}.}
#' Theorem 6 then sums the diagonal variances weighted by \eqn{a_i^2}
#' and twice the covariances weighted by \eqn{a_i a_l}. Every unordered
#' pair of entries is included, interpreting the theorem as the variance
#' of the weighted sum; the printed simultaneous coordinate-wise pair
#' bounds would omit crossed pairs. Powers and derivatives refer to the
#' restricted block, not the corresponding entries of the full \eqn{P^t}.
#' This reference uses a dense component covariance and is intended for
#' small operational blocks. The initial distribution is fixed; time zero
#' is extended to have variance zero.
#' @inherit vcov_rel return
#' @inherit rel references
#' @seealso [rel()], [vcov_rel()]
#' @inherit rel examples
#' @export
vcov_rel_sl <- function(x, time, up) {
  x <- validate_MarkovChain(x)
  validate_rel_time(time, scalar = TRUE)
  up <- validate_avail_up(up, get_state_names(x))
  Sigma <- rel_transition_covariance(x, up)

  r <- length(up)
  if (time == 0L || r == 0L || r == nstates(x)) return(0)

  powers <- Ppowers(x$P[up, up, drop = FALSE], max(1L, time - 1L))
  J <- matrix(0, r^2, r^2)
  for (i in seq_len(r)) {
    for (j in seq_len(r)) {
      for (u in seq_len(r)) {
        for (v in seq_len(r)) {
          for (h in seq_len(time)) {
            J[(i - 1L) * r + j, (u - 1L) * r + v] <-
              J[(i - 1L) * r + j, (u - 1L) * r + v] +
              powers[[h]][i, u] * powers[[time - h + 1L]][v, j]
          }
        }
      }
    }
  }
  S <- as.matrix(J %*% Sigma %*% t(J))
  weights <- rep(as.numeric(x$a[up]), each = r)
  variance <- sum(weights^2 * diag(S))
  for (j in seq_along(weights)) {
    if (j > 1L) {
      previous <- seq_len(j - 1L)
      variance <- variance + 2 * weights[j] *
        sum(weights[previous] * S[j, previous])
    }
  }
  as.numeric(variance)
}

# Select row-wise coordinates without constructing a dense full covariance.
rel_transition_covariance <- function(x, up) {
  positions <- match(up, get_state_names(x))
  coordinates <- as.vector(t(outer((positions - 1L) * nstates(x),
                                   positions, "+")))
  transition_covariance_sparse(x$P)[coordinates, coordinates, drop = FALSE]
}

validate_rel_time <- function(time, scalar = FALSE) {
  if (!is.numeric(time) || !is.null(dim(time)) || length(time) == 0L ||
      (scalar && length(time) != 1L) || anyNA(time) ||
      any(!is.finite(time)) || any(time < 0 | time != floor(time)) ||
      any(time >= .Machine$integer.max)) {
    requirement <- if (scalar) "a single" else "a nonempty vector of"
    stop(paste("`time` must be", requirement,
               "nonnegative integers less than .Machine$integer.max."),
         call. = FALSE)
  }
  invisible(NULL)
}
