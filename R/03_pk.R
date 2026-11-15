# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# k-step Transition Matrix                                                  ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Asymptotic Covariance Matrix for k-Step Transition Probabilities
#'
#' @description
#' Computes the asymptotic variance-covariance matrix of the row-wise
#' vectorized estimator of the \eqn{k}-step transition matrix \eqn{P^k} for an
#' irreducible finite-state Markov chain.
#'
#' @param x A `MarkovChain` or `StochMat` object representing an
#'   \eqn{s \times s} transition matrix \eqn{P}. The associated Markov chain
#'   must be irreducible.
#' @param k A positive integer specifying the transition step. Thus,
#'   \eqn{k = 1} corresponds to the one-step transition matrix, \eqn{k = 2}
#'   to the two-step transition matrix, and so forth.
#' @param ... Additional arguments passed to methods. Currently unused.
#'
#' @details
#' Let \eqn{\widehat P(n)} denote the empirical estimator of the one-step
#' transition matrix \eqn{P}. Then, for a fixed positive integer \eqn{k},
#' \deqn{
#'   \sqrt{n}
#'   \left\{
#'     \widehat P^k(n) - P^k
#'   \right\}
#'   \mathrel{\mathop{\longrightarrow}\limits^{d}}
#'   W_{P^k},
#' }
#'
#' For which,
#' \deqn{
#'   W_{\mathbf p^{(k)}}
#'   =
#'   \mathrm{vec}\left(W_{P^k}^{\top}\right)
#'   \sim
#'   \mathcal N_{s^2}
#'   \left(
#'     \mathbf 0_{s^2},
#'     \Sigma_{\mathbf p^{(k)}}
#'   \right).
#' }
#' This function calculates the asymptotic variance covariance matrix
#' \deqn{
#'   \Sigma_{\mathbf p^{(k)}}
#'   =
#'   M_k \Sigma_{\mathbf p} M_k^{\top},
#' }
#'
#' where
#' \deqn{
#'   M_k
#'   =
#'   \sum_{i=0}^{k-1}
#'   P^i \mathbin{\otimes}
#'   \left(P^{k-1-i}\right)^{\top},
#' }
#' and \eqn{\otimes} denotes the Kronecker product.
#'
#' The one-step covariance matrix \eqn{\Sigma_{\mathbf p}} is computed by
#' [vcov_p()].
#'
#' @return
#' A numeric matrix of dimension \eqn{s^2 \times s^2} containing
#' \eqn{\Sigma_{\mathbf p^{(k)}}}, where \eqn{s} is the number of states.
#'
#' @references
#' Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026).
#' \emph{Matrix asymptotic calculus for plug-in maximum likelihood estimators in finite Markov chains}.
#' arXiv preprint arXiv:2607.17161. \url{https://arxiv.org/abs/2607.17161}.
#' See Proposition 4 for the matrix differential and its Kronecker representation.
#'
#' Sadek, A., and Limnios, N. (2002).
#' \emph{Asymptotic properties for maximum likelihood estimators for reliability and failure rates of Markov chains}.
#' Communications in Statistics - Theory and Methods, 31(10), 1837--1861.
#' \doi{10.1081/STA-120014916}.
#' See Theorem 4 (p. 1847) for the asymptotic covariance matrix of the k-step
#' transition probabilities.
#'
#' @examples
#' P <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE))
#' dim(vcov_pk(P, k = 2L))
#' identical(vcov_pk(P, k = 1L), vcov_p(P))
#'
#' @seealso
#' [vcov_p()] for the one-step asymptotic covariance matrix.
#'
#' @export
vcov_pk <- function(x, k, ...) {
  UseMethod("vcov_pk", x)
}

#' @rdname vcov_pk
#' @export
vcov_pk.MarkovChain <- function(x, k, ...) {
  x <- validate_MarkovChain(x)
  method <- validated_stochmat_method(x$P, vcov_pk_validated, vcov_pk)
  method(x$P, k = k, ...)
}

#' @rdname vcov_pk
#' @importFrom fastmatrix kronecker.prod
#' @export
vcov_pk.StochMat <- function(x, k, ...) {
  vcov_pk_validated(validate_StochMat(x), k = k, ...)
}

vcov_pk_validated <- function(x, k, ...) {
  if (length(k) != 1L || !is.null(dim(k)) || !is_positive_integer(k)) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }
  k <- as.integer(k)

  if (k == 1L) {
    method <- validated_stochmat_method(x, vcov_p_validated, vcov_p)
    return(method(x))
  }

  Gamma <- transition_covariance_sparse(x)
  s <- nstates(x)
  Pk <- Ppowers(x, k - 1L) # P ^ 0, P ^ 1, ..., P ^ {k - 1}
  Mk <- matrix(0, s ^ 2, s ^ 2)
  for (i in seq_len(k)) {
    Mk <- Mk + fastmatrix::kronecker.prod(Pk[[i]], t(Pk[[k - i + 1L]]))
  }
  tcrossprod(as.matrix(Mk %*% Gamma), Mk)
}

#' @rdname vcov_pk
#' @export
vcov_pk.default <- function(x, k, ...) {
  stop("`x` must be a `MarkovChain` or `StochMat` object.", call. = FALSE)
}

#' Confidence Intervals for k-Step Transition Probabilities
#'
#' Compute a two-sided, pointwise normal-approximation confidence interval
#' for each entry of the k-step transition matrix \eqn{P^k}.
#'
#' @inheritParams ci_p
#' @param k A single positive integer giving the number of steps.
#'
#' @details
#' The estimates are the entries of the fitted transition matrix raised to
#' the kth power. The intervals use the diagonal of [vcov_pk()], based on
#' Proposition 4, Equation (26), of Gavrilopoulos, Trevezas, and Votsi (2026).
#' The normal approximation follows Remark 3, as in [ci_p()]. For `k = 1`,
#' the two functions give the same intervals.
#'
#' @inheritSection ci_p Pointwise coverage
#' @inherit ci_p return references
#' @seealso [ci_p()], [vcov_pk()], [Ppowers()]
#'
#' @examples
#' P <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE))
#' ci_pk(P, k = 3, n = 1000)
#' @export
ci_pk <- function(x, k, n, level = 0.95) {
  validate_ci_parameters(n, level)
  P <- ci_transition_matrix(x)
  variance <- diag(vcov_pk(x, k))
  Pk <- Ppowers(P, k)[[k + 1L]]

  # Structural zeros and rows with one reachable state have zero variance.
  # Enforce this before taking square roots to remove covariance roundoff.
  deterministic <- Pk == 0
  deterministic[rowSums(Pk > 0) == 1L, ] <- TRUE
  variance[as.vector(t(deterministic))] <- 0
  transition_ci(Pk, variance, n, level)
}

#' @title Reference Implementation of the k-Step Asymptotic Covariance Matrix
#'
#' @description
#' Computes the asymptotic variance-covariance matrix of the row-wise
#' vectorized estimator of the \eqn{k}-step transition matrix \eqn{P^k} for an
#' irreducible finite-state Markov chain.
#'
#' @param x A `MarkovChain` or `StochMat` object representing an
#'   \eqn{s \times s} transition matrix \eqn{P}. The associated Markov chain
#'   must be irreducible.
#' @param k A positive integer specifying the transition step. Thus,
#'   \eqn{k = 1} corresponds to the one-step transition matrix, \eqn{k = 2}
#'   to the two-step transition matrix, and so forth.
#'
#' @details
#' Unlike [vcov_pk()], which uses a matrix and Kronecker-product
#' implementation, this function follows the corresponding scalar or
#' component-wise expressions. It is used for internal testing only.
#'
#' @return A dense numeric matrix of dimension \eqn{s^2 \times s^2}
#'   containing the asymptotic variance-covariance matrix
#'   \eqn{\Sigma_{\mathbf p(k)}}, where \eqn{s} is the number of states.
#'
#' @references
#' Sadek, A., and Limnios, N. (2002).
#' \emph{Asymptotic properties for maximum likelihood estimators for reliability and failure rates of Markov chains}.
#' Communications in Statistics - Theory and Methods, 31(10), 1837--1861.
#' \doi{10.1081/STA-120014916}.
#' See Theorem 4 (p. 1847) for the component-wise covariance construction.
#'
#' @seealso [vcov_pk()] for the primary implementation and [vcov_p()] for the
#'   one-step asymptotic covariance matrix.
#'
#' @examples
#' P <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE))
#' all.equal(markov:::vcov_pk_sl(P, k = 2L), vcov_pk(P, k = 2L))
#'
#' @keywords internal
vcov_pk_sl <- function(x, k) {
  if (length(k) != 1L || !is.null(dim(k)) || !is_positive_integer(k)) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }
  k <- as.integer(k)
  if (inherits(x, "MarkovChain")) {
    x <- validate_MarkovChain(x)$P
  } else {
    if (!inherits(x, "StochMat")) {
      stop("`x` must be a `MarkovChain` or `StochMat` object.", call. = FALSE)
    }
    x <- validate_StochMat(x)
  }
  method <- validated_stochmat_method(x, vcov_p_validated, vcov_p)
  Gamma <- method(x)

  if (k == 1) {
    return(Gamma)
  } else {

    # Preliminaries
    s <- nstates(x)
    Mk <- matrix(0, nrow = s ^ 2, ncol = s ^ 2)
    Pk <- Ppowers(x, k - 1)[-1] # P ^ 1, ..., P ^ {k - 1}

    for (i in 1:s) {
      for (j in 1:s) {
        for (u in 1:s) {
          for (v in 1:s) {
            n <- (i - 1) * s + j
            l <- (u - 1) * s + v
            Mk[n, l] <- 1 * (j == v) * Pk[[k - 1]][i, u] +
              1 * (i == u) * Pk[[k - 1]][v, j]
            if (k > 2L) {
              a <- 0
              for (m in 2:(k - 1)) {
                a <- a + Pk[[k - m]][i, u] * Pk[[m - 1]][v, j]
              }
              Mk[n, l] <- Mk[n, l] + a
            }
          }
        }
      }
    }

    return(Mk %*% Gamma %*% t(Mk))
  }
}
