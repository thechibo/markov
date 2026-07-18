# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parallel Estimation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Parallel Estimation
#'
#' @description Simulate `k` trajectories of length `m` and apply `fun`.
#'
#' @param k integer. The number of trajectories.
#' @param m integer. The length of trajectories.
#' @param x an object of class `MarkovChain`.
#' @param fun the function to apply on each estimated MarkovChain.
#' @param ... extra arguments.
#'
#' @return A vector holding the combined return values of `fun` for each chain.
#'
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach getDoParWorkers foreach
#' @importFrom progress progress_bar
#' @export
#'
#' @inherit MarkovChain examples
estim_parallel <- function(x, fun, k = 1e3, m = 500, ...) {

  # Preliminaries
  traj <- list()
  y <- list()
  estim <- list()

  # Create a cluster
  n.cores <- parallel::detectCores() - 1
  my_cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  doSNOW::registerDoSNOW(cl = my_cluster)
  foreach::getDoParWorkers()

  # Progress Bar
  pb <- progress::progress_bar$new(
          format = "Processing [:bar] :elapsed | eta: :eta",
          total = k,
          width = 60)
  progress <- function() {
    pb$tick()
  }
  opts <- list(progress = progress)

  # Estimation for k trajectories
  estim <- foreach::foreach(i = 1:k, .combine = 'c', .options.snow = opts) %dopar% {
    traj[[i]] <- markov::sim(m, x)
    y[[i]] <- markov::fit(traj[[i]])
    y[[i]]$a <- x$a
    do.call(fun, list(x = y[[i]], ...))
  }

  # Stop the cluster
  parallel::stopCluster(my_cluster)

  # Return the estimation
  estim

}
