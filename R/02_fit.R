# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Estimate Markov Chain Parameters                                          ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Estimate a Markov Chain from One or More Trajectories
#'
#' @description
#' Estimates a time-homogeneous, discrete-time Markov chain from an observed
#' sequence of states or independent trajectories sharing a transition matrix.
#'
#' @param x An atomic vector or factor of observed states in time order, or
#'   an atomic matrix with one independent trajectory per row and observation
#'   times in columns. Each trajectory must have at least two observations;
#'   matrices must have at least one row. Missing values are not allowed.
#'   State labels must be non-empty.
#' @param init Initial distribution to attach to the fitted Markov chain.
#'   Must be exactly one of:
#'   \describe{
#'     \item{`"uniform"`}{Uniform distribution over the observed state space.}
#'     \item{`"realized"`}{Empirical distribution of the initial states across paths; for one path, concentrated on its first state. This is the default choice.}
#'     \item{`"stationary"`}{Stationary distribution implied by the estimated transition matrix (via `stat()`).}
#'   }
#'
#' @details
#' The state space consists of the distinct observed values, reading matrix
#' rows in order. If `x` is numeric or ordered, the
#' states are sorted; otherwise, the first-appearance order is used.
#' Numeric state labels use additional precision when needed to distinguish
#' values with identical default character representations.
#'
#' The function computes the empirical transition frequencies and normalizes by
#' row to obtain a transition probability matrix. Counts are pooled over paths;
#' no transitions are counted between the end of one path and the next path.
#'
#' If a state never appears as an origin (i.e., it is never observed at times
#' `t = 1, ..., n-1`), the corresponding transition row has no empirical counts.
#' In that case the function treats that state as absorbing (probability 1 on the
#' diagonal) and emits a warning.
#'
#' The row-normalized transition matrix is validated by [StochMat()].
#'
#' @return A `MarkovChain` object.
#'
#' @seealso [MarkovChain()], [sim()], [stat()]
#'
#' @inherit MarkovChain examples
#' @export
fit <- function(x, init = "realized") {

  if (!is.character(init) || length(init) != 1L ||
      !is.null(dim(init)) || !init %in% c("realized", "uniform", "stationary")) {
    stop(
      '`init` must be one of "realized", "uniform", or "stationary".',
      call. = FALSE
    )
  }

  x <- validate_trajectories(x)
  paths <- trajectory_data(x)

  # States (keep stable ordering for numeric/ordered)
  states <- unique(paths$values)
  if (is.numeric(states) || is.ordered(states)) {
    states <- sort(states)
  }
  s <- length(states)

  # Match values before formatting labels, preserving distinct numeric states.
  labels <- as.character(states)
  if (is.numeric(states) && anyDuplicated(labels)) {
    labels <- sprintf("%.17g", states)
  }
  if (anyNA(labels) || any(labels == "") || anyDuplicated(labels)) {
    stop("State labels in `x` must be unique, non-missing, and non-empty.",
         call. = FALSE)
  }
  xf <- factor(match(paths$values, states), levels = seq_len(s), labels = labels)

  # Transition counts
  from <- xf[paths$from]
  to   <- xf[paths$to]

  P <- unclass(table(from, to, dnn = NULL))

  # Row-normalize to transition probabilities
  rs <- rowSums(P)

  # Handle states that never appear as an origin (row sum = 0)
  if (any(rs == 0)) {
    zero_states <- labels[rs == 0]
    warning(
      "Some states never appear as an origin in `x`;
      treating them as absorbing states: ",
      paste(zero_states, collapse = ", ")
    )
    for (i in which(rs == 0)) {
      P[i, i] <- 1
    }
    rs <- rowSums(P)
  }

  P <- sweep(P, 1, rs, "/")

  # Initial distribution
  if (init == "uniform") {
    a <- rep(1 / s, s)
    names(a) <- labels
  } else if (init == "realized") {
    a <- setNames(tabulate(as.integer(xf[paths$first]), nbins = s) /
                    length(paths$first), labels)
  } else { # init == "stationary"
    a <- stat(StochMat(P))
    a <- setNames(as.numeric(a), labels)
  }

  MarkovChain(a = a, P = P)
}
