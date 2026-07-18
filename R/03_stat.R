# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stationary Distribution                                                   ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Stationary Distribution of a Markov Chain
#'
#' @description
#' Computes stationary distribution(s) of a finite-state, time-homogeneous
#' Markov chain.
#'
#' @param x A `MarkovChain` or a `StochMat` object.
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' For a transition matrix \eqn{P}, a stationary distribution \eqn{\pi}
#' satisfies \deqn{\pi^\top P = \pi^\top,\qquad \sum_i \pi_i = 1,\ \pi_i \ge 0.}
#'
#' If the chain has a single closed communicating class, the stationary
#' distribution is unique (though the chain may be periodic). If there are
#' multiple closed communicating classes, stationary distributions are not
#' unique. Every stationary distribution is a convex combination of the extremal
#' ones, each supported on one closed class.
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
  stat(x$P, ...)
}

#' @export
stat.StochMat <- function(x, ...) {

  # Preliminaries
  s <- nrow(x)
  st_names <- get_state_names(x)

  # Find extremals (stationary distributions for each recurrent class)
  closed_classes <- Filter(function(c) {is_rec_class(x, c)}, mc_classes(x))
  extremals <- lapply(closed_classes, function(c) {stat_on_class(x, c)})

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
#' The implementation uses the delta method in matrix form and Kronecker:
#' \deqn{\Sigma_\pi = M\,\Gamma\,M^\top,}
#' where \eqn{\Gamma} is the asymptotic covariance matrix of \eqn{\mathrm{vec}(\hat P)}
#' returned by [vcov_p()], and
#' \deqn{M = \pi^\top \otimes (I - P^\top + \mathbf{1}\mathbf{1}^\top)^{-1}.}
#'
#' @return A numeric \eqn{s \times s} matrix, where \eqn{s} is the number of
#' states, giving the asymptotic covariance matrix \eqn{\Sigma_\pi}.
#'
#' @references
#' GSI
#'
#' @seealso [stat()], [vcov_p()], [is_irreducible()]
#'
#' @inherit MarkovChain examples
#'
#' @export
vcov_stat <- function(x, ...) {
  UseMethod("vcov_stat", x)
}

#' @rdname vcov_stat
#' @export
vcov_stat.MarkovChain <- function(x, ...) {
  vcov_stat(x$P, ...)
}

#' @rdname vcov_stat
#' @export
vcov_stat.StochMat <- function(x, ...) {

  # Preliminaries
  if (!is_irreducible(x)) {
    stop("The Markov Chain must be irreducible.")
  }
  s <- nstates(x)
  Is <- diag(s)
  A <- matrix(1, nrow = s, ncol = s)
  pi_vec <- stat(x)

  # Delta-method transformation matrix
  M <- fastmatrix::kronecker.prod(matrix(pi_vec, nrow = 1L),
                                  solve(Is - t(x) + A))

  tcrossprod(M %*% vcov_p(x), M)
}

#' @rdname vcov_stat
#' @export
vcov_stat.default <- function(x, ...) {
  stop(
    "`x` must be a `MarkovChain` or `StochMat` object.",
    call. = FALSE
  )
}
