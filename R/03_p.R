# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Transition Probability Matrix                                             ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Row-Wise Asymptotic Covariance Block of Transition Probabilities
#'
#' @description
#' Computes the \eqn{m \times m} covariance block
#' \deqn{
#'   \Lambda_i = \operatorname{diag}(p_{i\bullet}) -
#'   p_{i\bullet}p_{i\bullet}^{\mathsf T},
#' }
#' associated with row \eqn{i} of a finite-state Markov chain transition
#' matrix \eqn{P}. Here,
#' \eqn{p_{i\bullet} = (p_{i1}, \ldots, p_{im})^{\mathsf T}} is the vector of
#' transition probabilities from state \eqn{i}.
#'
#' @param p Numeric probability row vector \eqn{p_{i\bullet}}. Should be
#'   non-negative and sum to 1 (row of a row-stochastic transition matrix).
#'
#' @details
#' This matrix is the covariance matrix of a single multinomial trial with
#' category-probability vector \eqn{p_{i\bullet}}. It also forms the
#' \eqn{i}-th diagonal block of the asymptotic covariance matrix of the
#' empirical transition-probability MLE \eqn{\widehat{\mathbf P}(n)}.
#'
#' The \eqn{(j,k)} entry of the matrix is
#' \deqn{
#'   \Lambda_i(j,k)
#'   =
#'   p_{ij}\left(\delta_{jk} - p_{ik}\right),
#'   \qquad j,k = 1,\ldots,m,
#' }
#' where \eqn{\delta_{jk}} is the Kronecker delta.
#'
#' @return A numeric matrix of dimension \eqn{s \times s}, where
#' \eqn{s = \mathrm{length}(p)} is the number of states.
#'
#' @seealso [vcov_p()]
#'
#' @keywords internal
lambda_block <- function(p) {
  p <- as.numeric(p)
  diag(p) - tcrossprod(p) # Lambda_i = diag(p) - p p^T
}

#' @title Asymptotic Variance–Covariance of the Transition Probability Vector
#'
#' @description
#' Constructs the block-diagonal asymptotic variance–covariance matrix for the
#' row-wise vectorization of an empirical transition matrix from an
#' irreducible finite-state Markov chain.
#'
#' @param x A `MarkovChain` or `StochMat` object representing an
#'   \eqn{s \times s} row-stochastic transition matrix. The Markov chain must
#'   be irreducible.
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' Let \eqn{\widehat{\mathbf P}(n)} be the empirical transition matrix
#' estimated from a trajectory of length \eqn{n}, and let
#' \eqn{\widehat{\mathbf p}(n)} be its row-wise vectorization. Under
#' irreducibility, the empirical transition matrix is strongly consistent and
#' asymptotically normal:
#' \deqn{
#'   \sqrt{n}
#'   \left\{
#'     \widehat{\mathbf p}(n) - \mathbf p
#'   \right\}
#'   \xrightarrow[n\to\infty]{d}
#'   \mathbf W_{\mathbf p},
#'   \qquad
#'   \mathbf W_{\mathbf p}
#'   \sim
#'   \mathcal N_{s^2}
#'   \left(
#'     \mathbf 0_{s^2},
#'     \Sigma_{\mathbf p}
#'   \right).
#' }
#'
#' This function calculates the matrix \eqn{\Sigma_{\mathbf p}}. Specifically,
#' if \eqn{\boldsymbol{\pi} = (\pi_1,\ldots,\pi_s)^{\mathsf T}} is the unique
#' stationary distribution, then the asymptotic covariance matrix is
#' \deqn{
#'   \Sigma_{\mathbf p}
#'   =
#'   \begin{pmatrix}
#'     \pi_1^{-1}\Lambda_1 & 0                         & \cdots & 0 \\
#'     0                   & \pi_2^{-1}\Lambda_2       & \cdots & 0 \\
#'     \vdots              & \vdots                    & \ddots & \vdots \\
#'     0                   & 0                         & \cdots &
#'       \pi_s^{-1}\Lambda_s
#'   \end{pmatrix},
#' }
#'
#' where
#' \deqn{
#'   \Lambda_i
#'   =
#'   \operatorname{diag}(p_{i\bullet})
#'   -
#'   p_{i\bullet}p_{i\bullet}^{\mathsf T}, \qquad i = 1,\ldots,s.
#' }
#'
#' @return
#' A numeric block-diagonal matrix of dimension \eqn{s^2 \times s^2},
#' where \eqn{s} is the number of states.
#'
#' @references
#' GSI - Theorem 1
#'
#' A. Sadek and N. Limnios (2002).
#' \emph{Asymptotic properties of estimators for nonhomogeneous Markov chains}.
#' Statistics & Probability Letters, 60(4), 1839–1850.
#' See Proposition 1 (p. 1846) for the asymptotic covariance matrix of the
#' transition probabilities.
#'
#' @inherit MarkovChain examples
#'
#' @importFrom Matrix bdiag
#' @export
vcov_p <- function(x, ...) {
  UseMethod("vcov_p", x)
}

#' @rdname vcov_p
#' @export
vcov_p.MarkovChain <- function(x, ...) {
  vcov_p(x$P, ...)
}

#' @rdname vcov_p
#' @export
vcov_p.StochMat <- function(x, ...) {
  if (!is_irreducible(x)) {
    stop("The Markov Chain must be irreducible.")
  }
  stat <- stat(x)

  blocks <- lapply(seq_len(nstates(x)), function(i) {
    lambda_block(x[i, ]) / stat[i]
  })

  as.matrix(Matrix::bdiag(blocks))
}


#' @rdname vcov_p
#' @export
vcov_p.default <- function(x, ...) {
  stop(
    "`x` must be a `MarkovChain` or `StochMat` object.",
    call. = FALSE
  )
}
