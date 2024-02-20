# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Classes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stochastic Vector
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

new_StochVec <- function(...) {
  x <- c(...)
  if (is.StochVec(x)) {
    return(x)
  }

  structure(x, class = c("StochVec", "numeric"))
}

validate_StochVec <- function(x) {

  if (any(is.na(x)) || any(x < 0)) {
    stop("All values must be non-missing and greater than zero.")
  }

  if (sum(x) != 1) {
    stop("Values must sum to 1.")
  }

  x

}

StochVec <- function(...) {
  validate_StochVec(new_StochVec(...))
}

is.StochVec <- function(x) {
  validate_StochVec(x)
  is(x, "StochVec")
}

setOldClass("StochVec")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stochastic Matrix
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

new_StochMat <- function(x, ...) {
  if (!is.matrix(x)) {
    x <- matrix(x, ...)
  }
  if (is.StochMat(x)) {
    return(x)
  }
  structure(x, class = c("StochMat", "matrix"))
}

validate_StochMat <- function(x) {

  if (!is.matrix(x) || (nrow(x) != ncol(x))) {
    stop("Input must be a square matrix.")
  }

  if (any(is.na(x)) || any(x < 0)) {
    stop("All values must be non-missing and greater than zero.")
  }

  if (!all(rowSums(x) == 1)) {
    stop("Row values must sum to 1.")
  }

  x

}

StochMat <- function(...) {
  validate_StochMat(new_StochMat(...))
}

is.StochMat <- function(x) {
  validate_StochMat(x)
  is(x, "StochMat")
}

setOldClass("StochMat")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Markov Chain
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

new_MarkovChain <- function(a, P) {
  x <- list(a = StochVec(a), P = StochMat(P))
  structure(x, class = c("MarkovChain", "list"))
}

validate_MarkovChain <- function(x) {

  if (!all(names(x) %in% c("a", "P"))) {
    stop("A Markov Chain must have two elements named 'a' and 'P'.")
  }

  is.StochVec(x$a)
  is.StochMat(x$P)

  x

}

#' MarkovChain S3 Class
#'
#' @param a numeric. The initial distribution vector.
#' @param P matrix. The transition probability matrix.
#'
#' @return An object of class 'MarkovChain' with two elements, named a and P.
#' @export
#'
#' @examples
#' # Choose the transition probability matrix and the initial law of the MC
#' a <- c(0.5, 0.5, 0)
#' P <- matrix(c(0.3, 0.7, 0,
#'               0.6, 0.1, 0.3,
#'               0.5, 0,   0.5), nrow = 3, ncol = 3, byrow = TRUE)
#'
#' # Create a Markov Chain
#' x <- MarkovChain(a = a, P = P)
#'
#' # Find Stationary distribution
#' pi_vec <- stationary(x)
#'
#' # Simulate a single trajectory of length m
#' m <- 1e4
#' traj <- sim(m, x)
#'
#' # Estimate the transition probability matrix
#' y <- estim(traj)
MarkovChain <- function(a, P) {
  validate_MarkovChain(new_MarkovChain(a, P))
}

is.MarkovChain <- function(x) {
  validate_MarkovChain(x)
  is(x, "MarkovChain")
}

setOldClass("MarkovChain")

new_MarkovChainList <- function(x) {
  stopifnot(is.list(x))
  structure(x, class = c("MarkovChainList", "list"))
}

validate_MarkovChainList <- function(x) {

  for (i in seq_along(x)) {
    if (!(is(x[[i]], "MarkovChain"))) {
      stop("All list elements be of class MarkovChain. Stopped at index", i)
    }
  }

  x

}

MarkovChainList <- function(x) {
  validate_MarkovChainList(new_MarkovChainList(x))
}

`[.MarkovChainList` <- function(x, i) {
  new_MarkovChainList(NextMethod())
}

setOldClass(c("MarkovChainList", "list"))
