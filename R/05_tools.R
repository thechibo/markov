# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tools
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# State Names                      ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @export
set_state_names <- function(x, ...) {
  UseMethod("set_state_names")
}

#' @export
set_state_names.MarkovChain <- function(x, ...) {
  x$a <- set_state_names(x$a, ...)
  x$P <- set_state_names(x$P, ...)
  x
}

#' @export
set_state_names.matrix <- function(x, ...) {
  if (is.null(rownames(x))) {
    rownames(x) <- as.character(seq_len(nrow(x)))
    colnames(x) <- as.character(seq_len(ncol(x)))
  }
  x
}

#' @export
set_state_names.numeric <- function(x, ...) {
  if (is.null(names(x))) {
    names(x) <- as.character(seq_along(x))
  }
  x
}

#' @export
get_state_names <- function(x, ...) {
  UseMethod("get_state_names")
}

#' @export
get_state_names.MarkovChain <- function(x, ...) {
  rownames(x$P)
}

#' @export
get_state_names.matrix <- function(x, ...) {
  rownames(x)
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
