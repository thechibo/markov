# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Hypothesis Tests                                                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Test a Specified Markov Transition Probability
#'
#' @description
#' Test the null hypothesis that a single transition probability equals
#' `p0`, against a two-sided alternative, using the asymptotic test of
#' Anderson and Goodman (1957).
#'
#' @param x An atomic vector or factor of states in time order, or an atomic
#'   matrix with one independent trajectory per row and common observation
#'   times in columns. Each path must have at least two observations; matrices
#'   must have at least one row. Missing values are not allowed. Validated by
#'   the internal `validate_trajectories()` helper.
#' @param ij An atomic vector or factor of length two, with no dimensions or
#'   missing values. Its first entry is the source state and its second is
#'   the destination state. Both must occur in `x`. Entries are state values,
#'   not positions in the state space; equal entries test a self-transition.
#' @param p0 A single finite numeric probability strictly between zero and
#'   one. Boundary nulls are excluded because the asymptotic null variance
#'   vanishes there.
#'
#' @details
#' Section 3.1 (pp. 96--97) of Anderson and Goodman (1957) derives tests for
#' specified transition probabilities from their multinomial asymptotic law.
#' For a single entry, aggregate destinations into state \eqn{j} and its
#' complement. Writing \eqn{N_i} for the number of departures from state
#' \eqn{i}, \eqn{N_{ij}} for transitions from \eqn{i} to \eqn{j}, and
#' \eqn{\widehat p_{ij}=N_{ij}/N_i}, the two-category version of Equation
#' (3.1) is
#' \deqn{X^2 = \frac{N_i(\widehat p_{ij}-p_0)^2}{p_0(1-p_0)}.}
#' Under the null, this has an asymptotic chi-squared distribution with one
#' degree of freedom. No continuity correction is applied. The final
#' observation is not counted as a departure; an error is raised if the
#' source state has no observed departures.
#'
#' The model is a time-homogeneous, first-order finite Markov chain. For a
#' single long trajectory, the approximation requires adequate visits to
#' the source state, as under the ergodic setting discussed in Section 5.
#' It is not an exact binomial test conditional on the random departure count.
#' For a matrix, rows are independent paths with a common transition matrix;
#' counts are pooled across paths, excluding each path's final observation
#' from departures. The Section 3.1 approximation applies as the number of
#' independent paths increases with their common length fixed.
#' A warning is issued if either null expected count, \eqn{N_i p_0} or
#' \eqn{N_i(1-p_0)}, is below five. The result reports a p-value rather than
#' imposing a significance level; reject at level \eqn{\alpha} when the
#' p-value is below \eqn{\alpha}.
#'
#' @return An object of class `htest`, printed using the standard R test
#'   format, with components `statistic` (named `X-squared`), `parameter`
#'   (one degree of freedom), `p.value`, `estimate`, `null.value`,
#'   `alternative` (`"two.sided"`), `method`, and `data.name`. Additional
#'   components `observed` and `expected` contain the two-category observed
#'   and null expected transition counts.
#'
#' @references
#' Anderson, T. W., and Goodman, L. A. (1957).
#' Statistical Inference about Markov Chains.
#' \emph{The Annals of Mathematical Statistics}, 28(1), 89--110.
#' \doi{10.1214/aoms/1177707039}.
#' See Section 3.1 for specified-probability tests and Section 5 for a
#' single long chain.
#' @seealso [fit()], [stats::chisq.test()], [stats::prop.test()]
#' @examples
#' x <- rep(c("a", "a", "b", "a", "b", "b"), 20)
#' test_prob(x, ij = c("a", "b"), p0 = 0.5)
#' @export
test_prob <- function(x, ij, p0) {
  data_name <- deparse1(substitute(x))
  x <- validate_trajectories(x)
  paths <- trajectory_data(x)
  if (!is.atomic(ij) || !is.null(dim(ij)) || length(ij) != 2L || anyNA(ij)) {
    stop("`ij` must contain exactly two non-missing state values.", call. = FALSE)
  }
  states <- unique(paths$values)
  pair <- match(ij, states)
  if (anyNA(pair)) {
    stop("Both values in `ij` must be states observed in `x`.", call. = FALSE)
  }
  if (!is.numeric(p0) || is.complex(p0) || !is.null(dim(p0)) || length(p0) != 1L ||
      !is.finite(p0) || p0 < 0 || p0 > 1) {
    stop("`p0` must be a single finite numeric probability in [0, 1].",
         call. = FALSE)
  }
  if (p0 == 0 || p0 == 1) {
    stop("The Anderson-Goodman test requires 0 < `p0` < 1.", call. = FALSE)
  }
  p0 <- as.numeric(p0)
  sequence <- match(paths$values, states)
  from <- sequence[paths$from] == pair[1L]
  departures <- sum(from)
  if (departures == 0) {
    stop("The source state in `ij` has no observed departures.", call. = FALSE)
  }
  successes <- sum(from & sequence[paths$to] == pair[2L])
  estimate <- successes / departures
  expected <- departures * c(p0, 1 - p0)
  if (any(expected < 5)) {
    warning("Chi-squared approximation may be incorrect: expected counts below 5.",
            call. = FALSE)
  }
  observed <- c(to_j = successes, other = departures - successes)
  statistic <- sum(((observed - expected) / sqrt(expected))^2)
  label <- paste0("p[", paste(as.character(ij), collapse = ", "), "]")
  structure(list(
    statistic = c("X-squared" = unname(statistic)),
    parameter = c(df = 1),
    p.value = stats::pchisq(statistic, df = 1, lower.tail = FALSE),
    estimate = setNames(estimate, label),
    null.value = setNames(as.numeric(p0), label),
    alternative = "two.sided",
    method = "Anderson-Goodman test of a transition probability",
    data.name = data_name,
    observed = observed,
    expected = setNames(as.numeric(expected), c("to_j", "other"))
  ), class = "htest")
}

#' Test Constancy of Markov Transition Probabilities across Time
#'
#' @description
#' Test whether independent Markov trajectories share a transition matrix
#' that is constant across observation times, using the Pearson test in
#' Section 3.2 of Anderson and Goodman (1957), Equation (3.8).
#'
#' @param x An atomic matrix with one independent trajectory per row and
#'   common observation times in columns. At least two rows and three columns
#'   are required. Missing values are not allowed. Validated by the internal
#'   `validate_trajectories()` helper. A single vector is not sufficient for
#'   this test.
#'
#' @details
#' Write \eqn{T = \mathrm{ncol}(x)-1} for the number of transition times and
#' \eqn{m} for the number of observed states. Let \eqn{N_{ij}(t)} count paths
#' making transition \eqn{i\to j} at time \eqn{t}, and let
#' \eqn{N_i(t-1)=\sum_j N_{ij}(t)}. Under the null, the pooled estimate is
#' \deqn{\widehat p_{ij} = \frac{\sum_t N_{ij}(t)}{\sum_t N_i(t-1)}.}
#' The statistic is
#' \deqn{X^2 = \sum_{i,j,t}
#'   \frac{(N_{ij}(t)-N_i(t-1)\widehat p_{ij})^2}
#'        {N_i(t-1)\widehat p_{ij}}.}
#' Its asymptotic null distribution is chi-squared with
#' \eqn{m(m-1)(T-1)} degrees of freedom. The alternative allows a different
#' transition matrix at each time. Paths must be independent, but observations
#' within a path need not be independent. The approximation is for increasing
#' numbers of paths with fixed \eqn{m} and \eqn{T}, positive transition
#' probabilities, and growing origin counts at every transition time.
#'
#' This implementation requires at least two observed states, a positive
#' origin count at every transition time for every state, and positive pooled
#' counts for every ordered state pair. Empty margins cause an error rather
#' than an automatic change to the paper's degrees of freedom. A warning is
#' issued when any expected cell count is below five. No continuity correction
#' is applied. Reject at a chosen significance level when the p-value is below
#' that level. This tests time homogeneity, not equality of the initial and
#' stationary distributions.
#'
#' @return An object of class `htest`, with components `statistic` (named
#'   `X-squared`), `parameter` (degrees of freedom), `p.value`, `method`, and
#'   `data.name`. Additional components are `estimate` (the pooled transition
#'   matrix), and `observed` and `expected` (arrays indexed by source state,
#'   destination state, and transition time). State order is first appearance
#'   when reading paths row by row.
#'
#' @references
#' Anderson, T. W., and Goodman, L. A. (1957).
#' Statistical Inference about Markov Chains.
#' \emph{The Annals of Mathematical Statistics}, 28(1), 89--110.
#' \doi{10.1214/aoms/1177707039}.
#' Section 3.2, especially Equations (3.6) and (3.8).
#' @seealso [test_prob()], [fit()], [stats::chisq.test()]
#' @examples
#' # Repeated copies of all binary paths of length three.
#' paths <- as.matrix(expand.grid(start = c("a", "b"),
#'                               middle = c("a", "b"), end = c("a", "b")))
#' x <- paths[rep(seq_len(nrow(paths)), each = 10), ]
#' test_const(x)
#' @export
test_const <- function(x) {
  data_name <- deparse1(substitute(x))
  x <- validate_trajectories(x)
  if (!is.matrix(x) || nrow(x) < 2L || ncol(x) < 3L) {
    stop("`x` must be a matrix with at least two paths and three observation times.",
         call. = FALSE)
  }
  states <- unique(as.vector(t(x)))
  m <- length(states)
  if (m < 2L) {
    stop("At least two observed states are required.", call. = FALSE)
  }
  times <- ncol(x) - 1L
  sequence <- matrix(match(x, states), nrow = nrow(x))
  labels <- as.character(states)
  if (is.numeric(states) && anyDuplicated(labels)) {
    labels <- sprintf("%.17g", states)
  }
  observed <- array(0, dim = c(m, m, times),
                    dimnames = list(from = labels, to = labels,
                                    time = as.character(seq_len(times))))
  for (t in seq_len(times)) {
    observed[, , t] <- matrix(tabulate(
      sequence[, t] + m * (sequence[, t + 1L] - 1L), nbins = m * m), m, m)
  }
  pooled <- apply(observed, c(1L, 2L), sum)
  origins <- apply(observed, c(1L, 3L), sum)
  if (any(origins == 0) || any(pooled == 0)) {
    stop(paste("The test requires positive origin counts at every time and",
               "positive pooled counts for every state pair."), call. = FALSE)
  }
  estimate <- pooled / rowSums(pooled)
  expected <- observed
  for (t in seq_len(times)) {
    expected[, , t] <- estimate * origins[, t]
  }
  if (any(expected < 5)) {
    warning("Chi-squared approximation may be incorrect: expected counts below 5.",
            call. = FALSE)
  }
  statistic <- sum(((observed - expected) / sqrt(expected))^2)
  df <- m * (m - 1) * (times - 1)
  structure(list(
    statistic = c("X-squared" = statistic),
    parameter = c(df = df),
    p.value = stats::pchisq(statistic, df = df, lower.tail = FALSE),
    method = "Anderson-Goodman test of constant transition probabilities",
    data.name = data_name,
    estimate = estimate,
    observed = observed,
    expected = expected
  ), class = "htest")
}
