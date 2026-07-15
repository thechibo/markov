# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Classes                                                                   ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stochastic Vector                ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname StochVec
#' @keywords internal
new_StochVec <- function(...) {
  x <- c(...)
  if (is.StochVec(x)) {
    return(x)
  }
  x <- set_state_names(x)
  structure(x, class = c("StochVec", "numeric"))
}

#' @rdname StochVec
#' @keywords internal
validate_StochVec <- function(x) {

  if (any(is.na(x)) || any(x < 0)) {
    stop("All values must be non-missing and greater than zero.")
  }

  if (sum(x) != 1) {
    stop("Values must sum to 1.")
  }

  x

}

#' @title Stochastic Vector Class
#'
#' @description
#' Construct and validate a stochastic vector: a numeric vector whose components
#' are strictly non-missing, strictly positive, and sum to exactly 1. The
#' returned object has S3 class `"StochVec"` and underlying storage `"numeric"`.
#'
#' @section Validation:
#' A valid `StochVec` must satisfy all of:
#' \itemize{
#'   \item no `NA`s
#'   \item all entries are `> 0`
#'   \item `sum(x) == 1`
#' }
#'
#' @param ... Numeric values coerceable to a flat numeric vector.
#'
#' @return
#' A `numeric` vector with additional class `"StochVec"`. If validation fails,
#' an error is thrown.
#'
#' @examples
#' StochVec(0.2, 0.3, 0.5)
#' x <- StochVec(c(0.1, 0.9))
#'
#' # check validity / class
#' is.StochVec(x)
#'
#' @seealso [new_StochVec()], [validate_StochVec()], [is.StochVec()]
#' @export
StochVec <- function(...) {
  validate_StochVec(new_StochVec(...))
}

#' @rdname StochVec
#' @export
is.StochVec <- function(x) {
  validate_StochVec(x)
  is(x, "StochVec")
}

setOldClass("StochVec")

#' @rdname StochVec
#' @export
print.StochVec <- function(x, digits = getOption("digits"), ...) {
  # ensure it's valid
  x <- validate_StochVec(x)

  cat("Stochastic Vector\n")
  cat("  length:", length(x), "\n")

  # drop class for numeric printing but keep names
  y <- unclass(x)
  y <- signif(y, digits = digits)
  if (!is.null(names(x))) {
    names(y) <- names(x)
  }

  # indent printed vector
  vec_lines <- capture.output(print(y, ...))
  indent <- "  "
  cat(paste0(indent, vec_lines), sep = "\n")

  invisible(x)
}

methods::setMethod("show", "StochVec", function(object) {
  print.StochVec(object)
})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stochastic Matrix                ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname StochMat
#' @keywords internal
new_StochMat <- function(x, tol = 1e-12, ...) {
  if (!is.matrix(x)) {
    x <- matrix(x, ...)
  }

  x <- validate_StochMat(x, tol = tol)
  x <- set_state_names(x)
  structure(x, class = c("StochMat", "matrix"))
}

#' @rdname StochMat
#' @keywords internal
validate_StochMat <- function(x, tol = 1e-12) {

  if (!is.matrix(x) || (nrow(x) != ncol(x))) {
    stop("Input must be a square matrix.")
  }
  if (!is.numeric(x)) {
    stop("Input must be a numeric matrix.")
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop("All values must be finite and non-missing.")
  }

  # Allow tiny negatives from floating point.
  if (any(x < - tol)) {
    stop("All values must be >= 0 (within numerical tolerance).")
  }
  if (any(x < 0)) {
    x[x < 0] <- 0
  }

  # Pure validation: accept if close enough, otherwise error.
  rs <- rowSums(x)
  close <- (abs(1 - rs) <= tol)
  if (!all(close)) {
    stop("Row values must sum to 1 (within numerical tolerance).")
  }

  x / rs
}

#' @title Stochastic Matrix Class
#'
#' @description
#' Construct and validate a stochastic matrix: a square numeric matrix whose
#' entries are strictly non-missing, strictly positive, and each row sums to
#' exactly 1. The returned object has S3 class `"StochMat"` and underlying
#' storage `"matrix"`.
#'
#' @section Validation:
#' A valid `StochMat` must satisfy all of:
#' \itemize{
#'   \item is a square matrix (`nrow(x) == ncol(x)`)
#'   \item no `NA`s
#'   \item all entries are `> 0`
#'   \item `rowSums(x) == 1` for every row
#' }
#'
#' @param ... Objects coercible to a matrix. Passed to `matrix()` when needed.
#' @param x For `new_StochMat()`, a matrix or something coercible to one.
#'
#' @return
#' A matrix with additional class `"StochMat"`. If validation fails,
#' an error is thrown.
#'
#' @examples
#' StochMat(diag(3))
#' StochMat(matrix(c(0.2,0.8, 0.7,0.3), ncol = 2, byrow = TRUE))
#'
#' m <- StochMat(diag(2))
#' is.StochMat(m)
#'
#' @seealso [new_StochMat()], [validate_StochMat()], [is.StochMat()]
#' @export
StochMat <- function(...) {
  validate_StochMat(new_StochMat(...))
}

#' @rdname StochMat
#' @export
is.StochMat <- function(x) {
  validate_StochMat(x)
  is(x, "StochMat")
}

setOldClass("StochMat")

#' @rdname StochMat
#' @export
print.StochMat <- function(x, digits = getOption("digits"), ...) {
  # ensure it's valid
  x <- validate_StochMat(x)

  cat("Stochastic Matrix\n")
  cat("  dim:", nrow(x), "x", ncol(x), "\n")

  # drop class for printing but keep dimnames
  y <- unclass(x)
  y[] <- signif(y, digits = digits)
  dimnames(y) <- dimnames(x)

  # indent printed matrix
  mat_lines <- capture.output(print(y, ...))
  indent <- "  "
  cat(paste0(indent, mat_lines), sep = "\n")

  invisible(x)
}

methods::setMethod("show", "StochMat", function(object) {
  print.StochMat(object)
})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Markov Chain                     ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname MarkovChain
#' @keywords internal
new_MarkovChain <- function(a, P) {
  x <- list(a = StochVec(a), P = StochMat(P))
  structure(x, class = c("MarkovChain", "list"))
}

#' @rdname MarkovChain
#' @keywords internal
validate_MarkovChain <- function(x) {
  if (!all(names(x) %in% c("a", "P"))) {
    stop("A Markov Chain must have two elements named 'a' and 'P'.")
  }
  is.StochVec(x$a)
  is.StochMat(x$P)
  x
}

#' Markov Chain (S3) Class
#'
#' Construct and validate a finite‐state, discrete‐time Markov chain, represented
#' by an initial distribution vector `a` and a row‐stochastic transition
#' probability matrix `P`. Both are validated and coerced to classes
#' [`StochVec`] and [`StochMat`].
#'
#' @section Structure:
#' A `MarkovChain` is a list with two named elements:
#' \describe{
#'   \item{a}{initial distribution, of class `StochVec`}
#'   \item{P}{transition matrix, of class `StochMat`}
#' }
#'
#' @section Validation:
#' A valid `MarkovChain` must satisfy:
#' \itemize{
#'   \item has exactly two elements named `"a"` and `"P"`
#'   \item `a` validates as `StochVec`
#'   \item `P` validates as `StochMat`
#' }
#'
#' @param a Numeric, initial law.
#' @param P Numeric matrix, transition probabilities (square row‐stochastic).
#'
#' @return
#' A list with class `"MarkovChain"` and elements `a` and `P`. Errors on
#' validation failure.
#'
#' @seealso [StochVec()], [StochMat()], [is.MarkovChain()]
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
#' pi_vec <- mc_stat(x)
#'
#' # Simulate a single trajectory of length m
#' m <- 1e4
#' traj <- mc_sim(m, x)
#'
#' # Estimate the transition probability matrix
#' y <- mc_estim(traj)
#'
#' # Make a list of Markov Chains
#' L <- MarkovChainList(list(x, y))
#' L[[1]]
MarkovChain <- function(a, P) {
  validate_MarkovChain(new_MarkovChain(a, P))
}

#' @rdname MarkovChain
#' @export
is.MarkovChain <- function(x) {
  validate_MarkovChain(x)
  is(x, "MarkovChain")
}

setOldClass("MarkovChain")

#' @rdname MarkovChain
#' @export
print.MarkovChain <- function(x, digits = getOption("digits"), ...) {
  # validate & unpack
  x <- validate_MarkovChain(x)
  a <- x$a
  P <- x$P
  n <- length(a)

  indent <- "    "  # 4 spaces

  cat("Markov Chain\n")
  cat("  Number of states:", n, "\n\n")

  ## Initial distribution
  cat("  Initial distribution (a):\n")
  a_num <- signif(unclass(a), digits = digits)
  if (!is.null(names(a))) names(a_num) <- names(a)

  a_lines <- capture.output(print(a_num, ...))
  cat(paste0(indent, a_lines), sep = "\n")
  cat("\n")

  ## Transition matrix
  cat("  Transition matrix (P):\n")
  P_mat <- unclass(P)
  P_mat[] <- signif(P_mat, digits = digits)
  dimnames(P_mat) <- dimnames(P)

  P_lines <- capture.output(print(P_mat, ...))
  cat(paste0(indent, P_lines), sep = "\n")

  invisible(x)
}

methods::setMethod("show", "MarkovChain", function(object) {
  print.MarkovChain(object)
})

#' @rdname MarkovChain
#' @keywords internal
new_MarkovChainList <- function(x) {
  stopifnot(is.list(x))
  structure(x, class = c("MarkovChainList", "list"))
}

#' @rdname MarkovChain
#' @keywords internal
validate_MarkovChainList <- function(x) {

  for (i in seq_along(x)) {
    if (!(is(x[[i]], "MarkovChain"))) {
      stop("All list elements be of class MarkovChain. Stopped at index", i)
    }
  }

  x

}

#' @rdname MarkovChain
#' @section MarkovChainList:
#' `MarkovChainList()` constructs a container of multiple `MarkovChain` objects.
#' All elements must be of class `"MarkovChain"`. Subsetting preserves the class.
#'
#' @param x For `MarkovChainList()`, a list of `MarkovChain` objects.
#'
#' @return For `MarkovChainList()`, an object of class `"MarkovChainList"`.
#'
#' @export
MarkovChainList <- function(x) {
  validate_MarkovChainList(new_MarkovChainList(x))
}

#' @rdname MarkovChain
#' @export
`[.MarkovChainList` <- function(x, i) {
  new_MarkovChainList(NextMethod())
}

setOldClass(c("MarkovChainList", "list"))
