# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tools
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# State Names                      ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Get or Set State Names
#'
#' @description
#' Utilities for assigning and retrieving the state names associated with
#' stochastic vectors, stochastic matrices, and Markov chains.
#'
#' `set_state_names()` assigns default sequential state names when names are
#' missing. Existing names are preserved.
#'
#' `get_state_names()` retrieves the state names attached to an object.
#'
#' @param x A `StochVec`, `StochMat`, or, for `get_state_names()`, a
#'   `MarkovChain` object.
#' @param ... Additional arguments passed to methods. Currently unused.
#'
#' @details
#' For an unnamed `StochVec`, `set_state_names()` assigns the names
#' `"1"`, `"2"`, ..., according to the vector length.
#'
#' For a `StochMat`, missing row or column names are assigned as
#' `"1"`, `"2"`, ..., according to the matrix dimensions. After assignment,
#' the row and column names must be identical. An error is raised when they
#' differ.
#'
#' For a `MarkovChain`, `get_state_names()` retrieves the row names of its
#' transition matrix.
#'
#' @return
#' `set_state_names()` returns `x` with state names assigned where necessary.
#'
#' `get_state_names()` returns a character vector containing the state names.
#'
#' @examples
#' a <- StochVec(c(0.4, 0.6))
#' get_state_names(a)
#'
#' P <- StochMat(matrix(c(0.8, 0.2,
#'                        0.3, 0.7),
#'                      nrow = 2, byrow = TRUE))
#' get_state_names(P)
#'
#' named_a <- StochVec(c(low = 0.4, high = 0.6))
#' get_state_names(named_a)
#'
#' named_P <- StochMat(
#'   matrix(c(0.8, 0.2,
#'            0.3, 0.7),
#'          nrow = 2, byrow = TRUE,
#'          dimnames = list(c("low", "high"), c("low", "high")))
#' )
#' get_state_names(named_P)
#'
#' mc <- MarkovChain(a = named_a, P = named_P)
#' get_state_names(mc)
#'
#' @name state_names
NULL

#' @rdname state_names
#' @export
set_state_names <- function(x, ...) {
  UseMethod("set_state_names", x)
}

#' @rdname state_names
#' @export
set_state_names.StochVec <- function(x, ...) {
  if (is.null(names(x))) {
    names(x) <- as.character(seq_along(x))
  }
  x
}

#' @rdname state_names
#' @export
set_state_names.StochMat <- function(x, ...) {
  if (is.null(rownames(x))) {
    rownames(x) <- as.character(seq_len(nrow(x)))
  }
  if (is.null(colnames(x))) {
    colnames(x) <- as.character(seq_len(ncol(x)))
  }
  if (!identical(rownames(x), colnames(x))) {
    stop("Row and column names must be identical.", call. = FALSE)
  }
  x
}

#' @rdname state_names
#' @export
set_state_names.default <- function(x, ...) {
  stop("`x` must be a `StochVec` or `StochMat` object.", call. = FALSE)
}

#' @rdname state_names
#' @export
get_state_names <- function(x, ...) {
  UseMethod("get_state_names")
}

#' @rdname state_names
#' @export
get_state_names.StochVec <- function(x, ...) {
  names(x)
}

#' @rdname state_names
#' @export
get_state_names.StochMat <- function(x, ...) {
  rownames(x)
}

#' @rdname state_names
#' @export
get_state_names.MarkovChain <- function(x, ...) {
  rownames(x$P)
}

#' @rdname state_names
#' @export
get_state_names.default <- function(x, ...) {
  stop("`x` must be a `StochVec`, a `StochMat`, or a `MarkovChain` object.",
       call. = FALSE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stationarity                     ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Partition states into communicating classes (equivalence classes)
#' @export
mc_classes <- function(x, ...) UseMethod("mc_classes")

#' @export
mc_classes.MarkovChain <- function(x, ...) {
  mc_classes(x$P, ...)
}

#' @export
mc_classes.StochMat <- function(x, ...) {
  s <- nrow(x)

  adj <- (x > 0)
  diag(adj) <- TRUE

  reach <- adj
  for (k in seq_len(s)) {
    reach <- reach | outer(reach[, k], reach[k, ], "&")
  }

  comm <- reach & t(reach)

  classes <- list()
  visited <- rep(FALSE, s)

  for (i in seq_len(s)) {
    if (visited[i]) next
    idx <- which(comm[i, ])
    visited[idx] <- TRUE
    classes[[length(classes) + 1L]] <- idx
  }

  classes
}

# Identify recurrent classes
is_recurrent <- function(P, c) {
  s <- nrow(P)
  outside <- setdiff(seq_len(s), c)
  if (length(outside) == 0L) {
    return(TRUE)
  }
  all(rowSums(P[c, outside, drop = FALSE]) <= 0)
}

# Compute extremal stationary distribution for each closed class
stat_on_class <- function(P, c) {

  P <- StochMat(P)
  s <- nrow(P)
  st_names <- get_state_names(P)
  Pc <- P[c, c, drop = FALSE]
  m <- length(c)

  # Solve (pi ^ T)(I - P) = 0 with sum(pi) = 1
  A <- rbind(t(Pc) - diag(1, m), rep(1, m))
  b <- c(rep(0, m), 1)

  pic <- as.numeric(qr.solve(A, b, tol = 1e-7))

  # Embed in full state space
  pi_full <- setNames(rep(0, s), st_names)
  pi_full[st_names[c]] <- pic / sum(pic)
  StochVec(pi_full)

}

#' Test Whether a Markov Chain Is Irreducible
#'
#' @description
#' Determines whether all states of a finite Markov chain communicate.
#'
#' @param x A `MarkovChain` or `StochMat` object.
#' @param ... Additional arguments passed to methods.
#'
#' @return A single logical value: `TRUE` if the chain is irreducible and
#'   `FALSE` otherwise.
#'
#' @export
is_irreducible <- function(x, ...) {
  UseMethod("is_irreducible", x)
}

#' @rdname is_irreducible
#' @export
is_irreducible.MarkovChain <- function(x, ...) {
  is_irreducible(x$P, ...)
}

#' @rdname is_irreducible
#' @export
is_irreducible.StochMat <- function(x, ...) {
  length(mc_classes(x)) == 1L
}

#' @rdname is_irreducible
#' @export
is_irreducible.default <- function(x, ...) {
  stop(
    "`x` must be a `MarkovChain` or `StochMat` object.",
    call. = FALSE
  )
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Kronecker                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cron <- function(i, n) {
  x <- rep(0, n)
  x[i] <- 1
  x
}

#' Number of States
#'
#' @description
#' Returns the number of states in a Markov chain or stochastic matrix.
#'
#' @param x A `MarkovChain` or `StochMat` object.
#'
#' @return A positive integer giving the number of states.
#'
#' @export
nstates <- function(x) {
  UseMethod("nstates", x)
}

#' @export
nstates.MarkovChain <- function(x) {
  nrow(x$P)
}

#' @export
nstates.StochMat <- function(x) {
  nrow(x)
}

#' @export
nstates.default <- function(x) {
  stop("`x` must be a `MarkovChain` or `StochMat` object.", call. = FALSE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check functions                  ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

is_positive_integer <- function(m) {
  valid <- length(m) == 1L &&
    is.numeric(m) &&
    !is.na(m) &&
    is.finite(m) &&
    m > 0 &&
    m == floor(m) &&
    m <= .Machine$integer.max

  if (!valid) {
    stop(sprintf("`%s` must be a positive integer.",
                 deparse(substitute(m))), call. = FALSE)
  }

  as.integer(m)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Matrix Calculus                  ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Matrix kronecker product
#' @export
kp <- function(x, y = x) {
  fastmatrix::kronecker.prod(x, y)
}

# Matrix cross product
#' @export
cp <- function(A, B) {
  tcrossprod(M %*% mat_gamma(x), M)
}


#' @export
Pdecomp <- function(P) {
  eigen(t(P))
}

#' @export
Ppower <- function (x, ...) {
  UseMethod("Ppower", x)
}

#' @export
Ppower.eigen <- function(x, k) {
  rvec <- x$vectors[, which.min(abs(x$values - 1))]
  rvec / sum(rvec)
  t(x$vectors %*% diag(x$values ^ k) %*% solve(x$vectors))
}

#' @export
Ppower.matrix <- function(x, k) {
  x <- Pdecomp(x)
  rvec <- x$vectors[, which.min(abs(x$values - 1))]
  rvec / sum(rvec)
  t(x$vectors %*% diag(x$values ^ k) %*% solve(x$vectors))
}

Ppowers <- function(P, k) {
  Pj <- vector("list", k)
  s <- nrow(P)
  Pj[[1L]] <- diag(1, s)              # P ^ 0
  for (j in 2L:(k + 1)) {
    Pj[[j]] <- Pj[[j - 1L]] %*% P   # P ^ (j - 1)
  }
  Pj
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parallel                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`%dopar%` <- foreach::`%dopar%`
`%do%` <- foreach::`%do%`
