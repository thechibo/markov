#' Simulate and Estimate Markov Chain Functionals
#'
#' Simulate independent trajectories, fit a Markov chain to each trajectory,
#' and apply a function to each fit.
#'
#' @param x A `MarkovChain` object.
#' @param fun A function accepting a fitted `MarkovChain` as argument `x`.
#' @param k A positive integer giving the number of trajectories.
#' @param m An integer of at least two giving the number of states observed per
#'   trajectory (`m - 1` transitions).
#' @param ... Additional arguments passed to `fun`.
#' @param n_workers A positive integer giving the maximum number of PSOCK
#'   workers. Set to one for sequential execution.
#'
#' @details
#' Each fit estimates the transition matrix with [fit()] and retains the
#' initial distribution from `x`. An error is raised if any trajectory omits
#' a state, because the fit would have a smaller state space. Choose a
#' trajectory length that gives adequate coverage of the states.
#'
#' Use [set.seed()] before calling this function for reproducible simulations
#' with a fixed worker count. For parallel execution, `fun` and its arguments
#' must be usable in a fresh R session with the package installed.
#'
#' @return The results of `fun` combined with [c()] in simulation order. For
#'   one trajectory, the result of `fun` is returned directly.
#' @seealso [sim()], [fit()]
#' @importFrom parallel clusterSetRNGStream makeCluster parLapply stopCluster
#' @export
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
#' x <- MarkovChain(c(1, 0), P)
#' set.seed(1)
#' estim_parallel(x, mttf, k = 4, m = 200, up = 1, n_workers = 1)
estim_parallel <- function(x, fun, k = 1e3, m = 500, ..., n_workers = 2L) {
  x <- validate_MarkovChain(x)
  fun <- match.fun(fun)
  check_integer <- function(value, name, minimum = 1L) {
    if (length(value) != 1L || !is.null(dim(value)) ||
        !is_positive_integer(value) || value < minimum) {
      stop(sprintf("`%s` must be a single integer >= %d.", name, minimum),
           call. = FALSE)
    }
  }
  check_integer(k, "k")
  check_integer(m, "m", minimum = 2L)
  check_integer(n_workers, "n_workers")

  iterations <- seq_len(k)
  n_workers <- min(n_workers, k)
  args <- list(...)
  if (n_workers == 1L) {
    results <- lapply(iterations, estimate_trajectory, chain = x, m = m,
                      callback = fun, args = args)
  } else {
    cluster <- parallel::makeCluster(n_workers, type = "PSOCK")
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterSetRNGStream(cluster, iseed = sample.int(.Machine$integer.max, 1L))
    results <- parallel::parLapply(cluster, iterations, estimate_trajectory,
                                   chain = x, m = m, callback = fun, args = args)
  }
  if (k == 1L) results[[1L]] else do.call(c, results)
}

estimate_trajectory <- function(i, chain, m, callback, args) {
  fitted <- fit(sim(m, chain))
  states <- get_state_names(chain)
  if (!all(states %in% get_state_names(fitted))) {
    stop(sprintf("Trajectory %d omitted at least one state.", i), call. = FALSE)
  }
  fitted <- MarkovChain(a = chain$a, P = fitted$P[states, states, drop = FALSE])
  do.call(callback, c(list(x = fitted), args))
}
