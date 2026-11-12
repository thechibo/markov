# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stationary Distribution                                                   ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Stationary Distribution of a Markov Chain
#'
#' @description
#' Computes stationary distribution(s) of a finite-state, time-homogeneous
#' Markov chain.
#'
#' @param x A valid `MarkovChain` or `StochMat` object.
#' @param ... Additional arguments passed to methods. Currently unused.
#'
#' @details
#' For a transition matrix \eqn{P}, a stationary distribution \eqn{\pi}
#' satisfies \deqn{\pi^\top P = \pi^\top,\qquad \sum_i \pi_i = 1,\ \pi_i \ge 0.}
#' For a `MarkovChain`, both components are validated and the stationary
#' distribution depends only on the transition matrix.
#'
#' If the chain has a single closed communicating class, the stationary
#' distribution is unique (though the chain may be periodic). If there are
#' multiple closed communicating classes, stationary distributions are not
#' unique. Every stationary distribution is a convex combination of the extremal
#' ones, each supported on one closed class.
#'
#' The Grassmann, Taksar and Heyman (1985) algorithm, also known as GTH,
#' computes the stationary probabilities within each closed class using
#' off-diagonal transitions.
#'
#' This function returns:
#' \itemize{
#'   \item a single `StochVec` if the stationary distribution is unique;
#'   \item otherwise a list of `StochVec` objects, containing extremal
#'   stationary distributions of the corresponding recurrent classes.
#' }
#'
#' @return
#' A `StochVec` (unique recurrent class), or a list of `StochVec` objects.
#'
#' @references
#' Grassmann, W. K., Taksar, M. I., and Heyman, D. P. (1985).
#' Regenerative Analysis and Steady State Distributions for Markov Chains.
#' \emph{Operations Research}, 33(5), 1107--1116.
#' \doi{10.1287/opre.33.5.1107}
#'
#' @seealso [MarkovChain()], [sim()], [fit()]
#'
#' @inherit MarkovChain examples
#'
#' @export
stat <- function(x, ...) {
  UseMethod("stat", x)
}

#' @export
stat.MarkovChain <- function(x, ...) {
  x <- validate_MarkovChain(x)
  method <- validated_stochmat_method(x$P, stat_validated, stat)
  method(x$P, ...)
}

#' @export
stat.StochMat <- function(x, ...) {
  stat_validated(validate_StochMat(x), ...)
}

stat_validated <- function(x, ...) {
  # Find extremals (stationary distributions for each recurrent class)
  closed_classes <- Filter(
    function(states) is_closed_class(x, states),
    mc_classes_validated(x)
  )
  extremals <- lapply(closed_classes, function(states) {
    stat_on_closed_class(x, states)
  })

  # Unique iff exactly one closed class
  if (length(closed_classes) == 1L) {
    return(extremals[[1L]])
  }

  extremals
}

#' @title Asymptotic Variance–Covariance Matrix of the stationary distribution MLE
#'
#' @description
#' Computes the asymptotic variance–covariance matrix of the maximum likelihood
#' estimator (MLE) of the stationary distribution of an irreducible,
#' finite-state Markov chain.
#'
#' @param x A `MarkovChain` or `StochMat` object. Must be irreducible.
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' Let \eqn{P} be the transition matrix of an irreducible finite-state Markov
#' chain with stationary distribution \eqn{\pi}. If \eqn{\hat P} is the MLE of
#' \eqn{P} based on a long trajectory, then the induced estimator
#' \eqn{\hat\pi = \pi(\hat P)} satisfies
#' \deqn{\sqrt{n}(\hat\pi - \pi) \xrightarrow{d} \mathcal{N}(0, \Sigma_\pi).}
#'
#' The implementation applies the matrix delta method using a Kronecker product:
#' \deqn{\Sigma_\pi = M\,\Gamma\,M^\top,}
#' where \eqn{\Gamma} is the asymptotic covariance matrix for the row-wise
#' vectorized transition estimator \eqn{\mathrm{vec}(\hat P^\top)} returned by
#' [vcov_p()]. With \eqn{\pi} written as a row vector,
#' \deqn{M = \pi \otimes (I - P^\top + \mathbf{1}\mathbf{1}^\top)^{-1}.}
#'
#' @return A numeric \eqn{s \times s} matrix, where \eqn{s} is the number of
#' states, giving the asymptotic covariance matrix \eqn{\Sigma_\pi}.
#'
#' @references
#' Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026).
#' \emph{Matrix asymptotic calculus for plug-in maximum likelihood estimators in finite Markov chains}.
#' arXiv preprint arXiv:2607.17161. \url{https://arxiv.org/abs/2607.17161}.
#' See Proposition 7 for the stationary differential and covariance.
#'
#' @seealso [stat()], [vcov_p()], [is_irreducible()]
#'
#' @examples
#' P <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE))
#' Sigma <- vcov_stat(P)
#' dim(Sigma) # One stationary probability per state
#'
#' @export
vcov_stat <- function(x, ...) {
  UseMethod("vcov_stat", x)
}

#' @rdname vcov_stat
#' @export
vcov_stat.MarkovChain <- function(x, ...) {
  x <- validate_MarkovChain(x)
  method <- validated_stochmat_method(x$P, vcov_stat_validated, vcov_stat)
  method(x$P, ...)
}

#' @rdname vcov_stat
#' @export
vcov_stat.StochMat <- function(x, ...) {
  vcov_stat_validated(validate_StochMat(x), ...)
}

vcov_stat_validated <- function(x, ...) {
  # Preliminaries
  pi_vec <- stationary_for_covariance(x)
  s <- nstates(x)
  Is <- diag(s)
  A <- matrix(1, nrow = s, ncol = s)

  # Delta-method transformation matrix
  M <- fastmatrix::kronecker.prod(matrix(pi_vec, nrow = 1L),
                                  solve(Is - t(x) + A))

  Gamma <- transition_covariance_sparse(x, stationary = pi_vec)
  tcrossprod(as.matrix(M %*% Gamma), M)
}

#' @rdname vcov_stat
#' @export
vcov_stat.default <- function(x, ...) {
  stop(
    "`x` must be a `MarkovChain` or `StochMat` object.",
    call. = FALSE
  )
}
