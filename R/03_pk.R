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
#' @param ... Additional arguments passed to [vcov_p()].
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
#'   \mathrm{vec}\left(W_{P^k}^{\mathMk T}\right)
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
#'   M_k \Sigma_{\mathbf p} M_k^{\mathMk T},
#' }
#'
#' where
#' \deqn{
#'   M_k
#'   =
#'   \sum_{i=0}^{k-1}
#'   P^i \mathbin{\otimes}
#'   \left(P^{k-1-i}\right)^{\mathMk T},
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
#' GSI
#'
#' A. Sadek and N. Limnios (2002).
#' \emph{Asymptotic properties of estimators for nonhomogeneous Markov chains}.
#' Statistics & Probability Letters, 60(4), 1839–1850.
#' See Theorem 4 (p. 1847) for the asymptotic covariance matrix of the k-step
#' transition probabilities.
#'
#' @inherit MarkovChain examples
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
  vcov_pk(x$P, k = k, ...)
}

#' @rdname vcov_pk
#' @importFrom fastmatrix kronecker.prod
#' @export
vcov_pk.StochMat <- function(x, k, ...) {

  if (k == 1) {
    return(vcov_p(x))
  } else {
    s <- nstates(x)
    Pk <- Ppowers(x, k - 1) # P ^ 0, P ^ 1, ..., P ^ {k - 1}
    Mk <- matrix(0, s ^ 2, s ^ 2)
    for (i in 1:k) {
      Mk <- Mk + fastmatrix::kronecker.prod(Pk[[i]], t(Pk[[k - i + 1]]))
    }
    return(Mk %*% vcov_p(x) %*% t(Mk))
  }
}

#' @rdname vcov_pk
#' @export
vcov_pk.default <- function(x, k, ...) {
  stop("`x` must be a `MarkovChain` or `StochMat` object.", call. = FALSE)
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
#' A. Sadek and N. Limnios (2002).
#' \emph{Asymptotic properties of estimators for nonhomogeneous Markov chains}.
#' Statistics & Probability Letters, 60(4), 1839--1850.
#'
#' @seealso [vcov_pk()] for the primary implementation and [vcov_p()] for the
#'   one-step asymptotic covariance matrix.
#'
#' @inherit MarkovChain examples
#'
#' @keywords internal
vcov_pk_sl <- function(x, k) {

  if (k == 1) {
    return(vcov_p(x))
  } else {

    # Preliminaries
    s <- nstates(x)
    Mk <- matrix(0, nrow = s ^ 2, ncol = s ^ 2)
    Pk <- Ppowers(x$P, k - 1)[-1] # P ^ 1, ..., P ^ {k - 1}

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

    return(Mk %*% vcov_p(x) %*% t(Mk))
  }
}
