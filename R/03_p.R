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
  block <- -tcrossprod(as.numeric(p))
  diag(block) <- 0
  # The diagonal is p[j] times the total probability of the other states.
  diag(block) <- -rowSums(block)
  block
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
#' Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026).
#' \emph{Matrix asymptotic calculus for plug-in maximum likelihood estimators in finite Markov chains}.
#' arXiv preprint arXiv:2607.17161. \url{https://arxiv.org/abs/2607.17161}.
#' See Theorem 1 for the row-wise transition covariance.
#'
#' Sadek, A., and Limnios, N. (2002).
#' \emph{Asymptotic properties for maximum likelihood estimators for reliability and failure rates of Markov chains}.
#' Communications in Statistics - Theory and Methods, 31(10), 1837--1861.
#' \doi{10.1081/STA-120014916}.
#' See Proposition 1 (p. 1846) for the asymptotic covariance matrix of the
#' transition probabilities.
#'
#' @examples
#' P <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE))
#' Sigma <- vcov_p(P)
#' dim(Sigma) # Four transition probabilities in row-wise order
#'
#' @importFrom Matrix bdiag
#' @export
vcov_p <- function(x, ...) {
  UseMethod("vcov_p", x)
}

#' @rdname vcov_p
#' @export
vcov_p.MarkovChain <- function(x, ...) {
  x <- validate_MarkovChain(x)
  method <- validated_stochmat_method(x$P, vcov_p_validated, vcov_p)
  method(x$P, ...)
}

#' @rdname vcov_p
#' @export
vcov_p.StochMat <- function(x, ...) {
  vcov_p_validated(validate_StochMat(x), ...)
}

vcov_p_validated <- function(x, ...) {
  as.matrix(transition_covariance_sparse(x))
}


#' @rdname vcov_p
#' @export
vcov_p.default <- function(x, ...) {
  stop(
    "`x` must be a `MarkovChain` or `StochMat` object.",
    call. = FALSE
  )
}

stationary_for_covariance <- function(x) {
  if (length(mc_classes_validated(x)) != 1L) {
    stop("The Markov Chain must be irreducible.", call. = FALSE)
  }
  stat_on_closed_class(x, get_state_names(x))
}

transition_covariance_sparse <- function(x, stationary = stationary_for_covariance(x)) {
  blocks <- lapply(seq_len(nstates(x)), function(i) {
    lambda_block(x[i, ]) / stationary[i]
  })
  Matrix::bdiag(blocks)
}

#' Confidence Intervals for Transition Probabilities
#'
#' Compute a two-sided, pointwise normal-approximation confidence interval
#' for each entry of the one-step transition matrix.
#'
#' @param x A [MarkovChain()] or [StochMat()] containing the transition
#'   matrix estimated from a single observed trajectory. The estimated
#'   chain must be irreducible.
#' @param n A single finite positive integer giving the number of observed
#'   one-step transitions used to estimate the transition matrix.
#' @inheritParams ci_rel
#'
#' @details
#' The intervals use the diagonal of the row-wise asymptotic covariance
#' matrix from [vcov_p()], as in Theorem 1, Equation (2), of Gavrilopoulos,
#' Trevezas, and Votsi (2026). The standard errors divide the square roots
#' of these variances by \eqn{\sqrt{n}} before applying the normal quantile
#' in Remark 3.
#'
#' @section Pointwise coverage:
#' The confidence level applies to each entry separately. Bounds may extend
#' outside \eqn{[0,1]}. A zero estimated variance gives a collapsed interval.
#' Nominal asymptotic coverage requires a positive limiting variance.
#'
#' @return A data frame with one row per transition, ordered by origin state
#'   and then by destination state in matrix order. Character columns
#'   `from` and `to` identify the states. Numeric columns `estimate`,
#'   `std.error`, `lower`, and `upper` give the plug-in probability,
#'   its estimated standard error, and its confidence limits.
#' @inherit ci_rel references
#' @seealso [vcov_p()], [ci_pk()]
#'
#' @examples
#' P <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE))
#' ci_p(P, n = 1000)
#' @export
ci_p <- function(x, n, level = 0.95) {
  validate_ci_parameters(n, level)
  P <- ci_transition_matrix(x)
  transition_ci(P, diag(vcov_p(x)), n, level)
}

ci_transition_matrix <- function(x) {
  if (inherits(x, "MarkovChain")) {
    return(validate_MarkovChain(x)$P)
  }
  if (inherits(x, "StochMat")) {
    return(validate_StochMat(x))
  }
  stop("`x` must be a `MarkovChain` or `StochMat` object.", call. = FALSE)
}

transition_ci <- function(P, variance, n, level) {
  intervals <- norm_ci(as.vector(t(P)), variance, n, level)
  data.frame(from = rep(rownames(P), each = ncol(P)),
             to = rep(colnames(P), times = nrow(P)), intervals)
}
