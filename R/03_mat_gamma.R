# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Gamma Matrix                                                              ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Row-wise asymptotic covariance block of transition probabilities
#'
#' @description
#' Compute the \eqn{m \times m} block
#' \deqn{\Gamma_i = \mathrm{diag}(p_{i\bullet}) - p_{i\bullet} p_{i\bullet}^\top}
#' corresponding to a single row \eqn{p_{i\bullet}} of a Markov chain transition
#' matrix \eqn{P}, where \eqn{p_{i\bullet}} denotes the row vector of
#' transition probabilities from state \eqn{i}.
#'
#' @param p Numeric probability row vector \eqn{p_{i\bullet}}. Should be
#'   non-negative and sum to 1 (row of a row-stochastic transition matrix).
#'
#' @return A numeric matrix of dimension \eqn{m \times m}, where
#'   \eqn{m = \mathrm{length}(p)}.
#'
#' @keywords internal
lambda_block <- function(p) {
  p <- as.numeric(p)
  diag(p) - tcrossprod(p) # Γ_i = diag(p) - p p^T
}

#' @title Asymptotic Variance–Covariance of the Transition Probability Vector
#'
#' @description
#' Constructs the block-diagonal asymptotic variance–covariance matrix of the
#' vectorized transition probability matrix of an irreducible finite-state
#' Markov chain.
#'
#' @param x A `MarkovChain` or `StochMat` object. It must be irreducible.
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' For row \eqn{i}, the block \eqn{\Gamma_i / \pi_i} is placed on the
#' diagonal, where
#' \eqn{\Gamma_i = \mathrm{diag}(p_{i\bullet}) - p_{i\bullet} p_{i\bullet}^\top}
#' and \eqn{\pi_i} is the stationary probability of state \eqn{i}.
#'
#' @return A numeric block-diagonal matrix of size \eqn{s^2 \times s^2},
#' where \eqn{s} is the number of states.
#'
#' @references
#' A. Sadek and N. Limnios (2002).
#' \emph{Asymptotic properties of estimators for nonhomogeneous Markov chains}.
#' Statistics & Probability Letters, 60(4), 1839–1850. (See p.1846 for Γ.)
#'
#' @inherit MarkovChain examples
#'
#' @export
mat_gamma <- function(x, ...) {
  UseMethod("mat_gamma", x)
}

#' @rdname mat_gamma
#' @export
mat_gamma.MarkovChain <- function(x, ...) {
  mat_gamma(x$P, ...)
}

#' @rdname mat_gamma
#' @export
mat_gamma.StochMat <- function(x, ...) {
  if (!is_irreducible(x)) {
    stop("The Markov Chain must be irreducible.")
  }
  stat <- mc_stat(x)

  blocks <- lapply(seq_len(nstates(x)), function(i) {
    lambda_block(x[i, ]) / stat[i]
  })

  as.matrix(Matrix::bdiag(blocks))
}


#' @rdname mat_gamma
#' @export
mat_gamma.default <- function(x, ...) {
  stop(
    "`x` must be a `MarkovChain` or `StochMat` object.",
    call. = FALSE
  )
}
