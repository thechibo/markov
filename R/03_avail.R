# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Availability
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Availability
#'
#' @param x an object of class `MarkovChain`.
#' @param time integer vector. The time-steps of interest.
#'
#' @return numeric. The probability of the system to be available at `time`.
#' @export
#'
#' @inherit MarkovChain examples
avail <- function(x, time) {

  nstate <- length(x$a)
  nup <- sum(x$a != 0)
  vec_ones <- rep(c(1, 0), c(nup, nstate - nup))

  y <- numeric(length(time))
  for (i in seq_along(time)) {
    y[i] <- as.numeric(x$a %*% Ppower(x$P, time[i]) %*% vec_ones)
  }

  y

}

#' @export
avail_gsi <- function(n, init, mylist, vec_ones, mat) {

  nstate <- nrow(mylist[[1]])
  s <- matrix(0, ncol = nstate ^ 2, nrow = 1)

  for (i in 1:n) {
    tmp <- fastmatrix::kronecker.prod(init %*% mylist[[i]],
                                      t(mylist[[n - i + 1]] %*% t(vec_ones)))
    s <- s + tmp
  }

  s %*% mat %*% t(s)

}

#This function returns the variance of the availability as a function
#of the time "n", Theorem 5, Sadek and Limnios
avail_var <- function(x, n, nup) {

  init <- rep(x$a, each = length(x$a))
  t(init) %*% avail_cov(x, n, nup) %*% init

}

avail_cov <- function(x, n, nup){

  mat_sigma <- sigma_mat(x, n)
  stat <- stationary(x)
  nstate <- length(stat)

  for (i in seq_along(stat)) {
    k <- (i - 1) * nstate + nup + 1
    mat_sigma[k, ] <- 0
    mat_sigma[, k] <- 0
  }

  mat_sigma

}
