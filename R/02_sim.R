# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate a Markov Chain Trajectory
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Markov Chain Trajectory Simulation
#'
#' @param m integer. The chain trajectory length.
#' @param x an object of class `MarkovChain`.
#'
#' @return A vector of length `m` with values simulated from the MarkovChain in
#' `x`.
#' @export
#'
#' @inherit MarkovChain examples
sim <- function(m, x) {

  # # Preliminaries
  # chain <- numeric(m)
  # states <- seq_along(x$a)
  #
  # # Simulate the initial state
  # chain[1] <- sample(states, size = 1, prob = x$a)
  #
  # for (i in 2:m) {
  #   chain[i] <- sample(states, size = 1, prob = x$P[chain[i - 1], ])
  # }
  #
  # chain


  ###########
  #
  chain          <- numeric(m)
  nstate         <- nrow(x$P)
  starting_value <- sample.int(nstate,size=1,replace=FALSE,prob=c(rep(1/nstate,nstate)))
  chain[1]       <- starting_value                      #simulate the initial state

  transition_probabilities <- x$P

  for (i in 2:m) {
    this_step_distribution  <- transition_probabilities[chain[i-1],]
    cumulative_distribution <- cumsum(this_step_distribution)
    r                      <- runif(1, min = 0, max = 1)
    chain[i]                <- min(which((cumulative_distribution>r) == TRUE))    #simulate the next state visited by the EMC
  }
  #
  chain
  #############



}
