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
  structure(x, class = c("StochVec", "numeric"))
}

#' @rdname StochVec
#' @keywords internal
validate_StochVec <- function(x, tol = 1e-12) {
  if (!is.numeric(x)) {
    stop("Input must be numeric.", call. = FALSE)
  }

  if (anyNA(x) || any(!is.finite(x))) {
    stop("All values must be finite and non-missing.", call. = FALSE)
  }

  if (any(x < 0)) {
    stop("All values must be nonnegative.", call. = FALSE)
  }

  total <- sum(x)
  if (abs(total - 1) > tol) {
    stop(
      sprintf("Values must sum to 1 within %g. Observed sum: %.17g.", tol, total),
      call. = FALSE
    )
  }

  # Set and check the state names
  x <- set_state_names(x)

  x / total
}

#' @title Stochastic Vector Class
#'
#' @description
#' Construct and validate a stochastic vector of finite, nonnegative
#' probabilities. The returned object has S3 class `"StochVec"`.
#'
#' @details
#' If the sum differs from one by at most `1e-12`, `StochVec()` divides
#' every entry by that sum without a warning. A larger deviation is an error.
#'
#' @section Validation:
#' A valid `StochVec` must satisfy all of:
#' \itemize{
#'   \item all entries are finite and nonnegative
#'   \item the sum is within `1e-12` of one
#'   \item state names are unique, non-missing, and non-empty
#' }
#'
#' @param ... For the constructors, numeric values combined into a flat vector.
#'   For the print method, additional arguments passed to [base::print()].
#' @param x An object to check with `is.StochVec()`, or a `StochVec` to validate
#'   or print.
#' @param tol Nonnegative numeric tolerance for the absolute difference between
#'   the sum and one. The validator defaults to `1e-12`.
#' @param digits Number of significant digits used to round displayed values.
#'   Defaults to `getOption("digits")`.
#'
#' @return
#' The constructors and validator return a numeric vector of class
#' `"StochVec"`. `is.StochVec()` returns a single logical value. The print method
#' invisibly returns the validated vector.
#'
#' @importFrom utils capture.output
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
  tryCatch(
    {
      validate_StochVec(x)
      inherits(x, "StochVec")
    },
    error = function(e) { FALSE }
  )
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
new_StochMat <- function(x, ...) {
  if (!is.matrix(x)) {
    x <- matrix(x, ...)
  }
  structure(x, class = c("StochMat", "matrix"))
}

#' @rdname StochMat
#' @keywords internal
validate_StochMat <- function(x, tol = 1e-12) {

  if (!is.matrix(x) || (nrow(x) != ncol(x))) {
    stop("Input must be a square matrix.", call. = FALSE)
  }
  if (nrow(x) == 0L) {
    stop("The state space must be non-empty.", call. = FALSE)
  }
  if (!is.numeric(x)) {
    stop("Input must be a numeric matrix.", call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop("All values must be finite and non-missing.", call. = FALSE)
  }

  # Allow tiny negatives from floating point.
  if (any(x < -1e-16)) {
    stop("All values must be >= 0 (within numerical tolerance).", call. = FALSE)
  }
  x[x < 0] <- 0

  # Check the row sums after clipping.
  rs <- rowSums(x)
  bad_row <- match(TRUE, abs(rs - 1) > tol, nomatch = 0L)
  if (bad_row != 0L) {
    stop(
      sprintf(
        "Row values must sum to 1 within %g. Row %d sums to %.17g.",
        tol, bad_row, rs[bad_row]
      ),
      call. = FALSE
    )
  }

  # Set and check the state names
  x <- set_state_names(x)

  x / rs
}

#' @title Stochastic Matrix Class
#'
#' @description
#' Construct and validate a square matrix of finite, nonnegative transition
#' probabilities. The returned object has S3 class `"StochMat"`.
#'
#' @details
#' Input entries below `-1e-16` are rejected. Entries in `[-1e-16, 0)`
#' are set to zero. After clipping, each row sum must differ from one by at most
#' `1e-12`. A larger deviation is an error. Each accepted row is divided by
#' its sum. `StochMat()` warns only if it clips negative entries.
#'
#' @section Validation:
#' A valid `StochMat` must satisfy all of:
#' \itemize{
#'   \item is a square matrix (`nrow(x) == ncol(x)`)
#'   \item has at least one state
#'   \item all entries are finite and nonnegative
#'   \item every row sum is within `1e-12` of one
#'   \item state names are unique, non-missing, and non-empty
#' }
#'
#' Assign the returned object to retain corrections. `is.StochMat()` only
#' checks validity.
#'
#' @param ... For the constructors, data and arguments for [base::matrix()]
#'   when the input is supplied as data rather than as a matrix. For the print
#'   method, additional arguments passed to [base::print()].
#' @param x For `new_StochMat()`, a matrix or data for [base::matrix()]. For
#'   `is.StochMat()`, the object to check. For validation or printing, a `StochMat`.
#' @param tol Nonnegative numeric tolerance for the absolute difference between
#'   each row sum and one after clipping. The validator defaults to `1e-12`.
#' @param digits Number of significant digits used to round displayed values.
#'   Defaults to `getOption("digits")`.
#'
#' @return
#' The constructors and validator return a matrix of class `"StochMat"`.
#' `is.StochMat()` returns a single logical value. The print method invisibly
#' returns the validated matrix.
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
  x <- new_StochMat(...)
  validated <- validate_StochMat(x)
  if (any(x < 0)) {
    warning("Tiny negative entries were set to zero and rows were normalized.",
            call. = FALSE)
  }
  validated
}

#' @rdname StochMat
#' @export
is.StochMat <- function(x) {
  tryCatch(
    {
      validate_StochMat(x)
      inherits(x, "StochMat")
    },
    error = function(e) { FALSE }
  )
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
  if (!inherits(x, "MarkovChain") || !is.list(x)) {
    stop("`x` must be a MarkovChain object.", call. = FALSE)
  }
  if (length(x) != 2L || !setequal(names(x), c("a", "P"))) {
    stop("A Markov Chain must have exactly two elements named 'a' and 'P'.",
         call. = FALSE)
  }
  a <- tryCatch(validate_StochVec(x$a), error = function(e) NULL)
  if (is.null(a) || !inherits(x$a, "StochVec")) {
    stop("`a` must be a valid StochVec object.", call. = FALSE)
  }
  P <- tryCatch(validate_StochMat(x$P), error = function(e) NULL)
  if (is.null(P) || !inherits(x$P, "StochMat")) {
    stop("`P` must be a valid StochMat object.", call. = FALSE)
  }
  if (length(x$a) != nrow(x$P)) {
    stop("The initial distribution and transition matrix must be of compatible dimensions.")
  }
  # Check the state names
  if (!identical(get_state_names(x$a), get_state_names(x$P))) {
    stop("The initial distribution and transition matrix must have identical state names.", call. = FALSE)
  }
  x$a <- a
  x$P <- P
  x
}

# Reuse validation for the package's StochMat class and dispatch on subclasses.
validated_stochmat_method <- function(P, implementation, generic) {
  if (identical(class(P), c("StochMat", "matrix"))) {
    return(implementation)
  }
  generic
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
#'   \item `length(a)` equals `nrow(P)`
#'   \item the names of `a` and the row and column names of `P` are identical,
#'     in the same order
#' }
#'
#' @param a Numeric vector of initial state probabilities, one per row of `P`.
#' @param P Square numeric row-stochastic matrix of transition probabilities.
#'   See [set_state_names()] for the default names assigned during construction.
#' @param x An object to check with `is.MarkovChain()`, a `MarkovChain` to
#'   validate or print, a list of valid chains to construct a `MarkovChainList`,
#'   or a `MarkovChainList` to subset.
#' @param digits Number of significant digits used to round displayed values.
#'   Defaults to `getOption("digits")`.
#' @param ... Additional arguments passed to [base::print()] for the initial
#'   distribution and transition matrix.
#'
#' @return
#' `MarkovChain()` returns a list of class `"MarkovChain"` with elements `a` and
#' `P`. `is.MarkovChain()` returns a single logical value. The print method
#' invisibly returns the validated chain.
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
#' pi_vec <- stat(x)
#'
#' # Simulate a single trajectory of length m
#' m <- 1e4
#' traj <- sim(m, x)
#'
#' # Estimate the transition probability matrix
#' y <- fit(traj)
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
  tryCatch(
    {
      validate_MarkovChain(x)
      inherits(x, "MarkovChain")
    },
    error = function(e) { FALSE }
  )
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
  if (!is.list(x)) {
    stop("`x` must be a list of valid MarkovChain objects.", call. = FALSE)
  }
  for (i in seq_along(x)) {
    if (!is.MarkovChain(x[[i]])) {
      stop("All list elements must be valid MarkovChain objects. Stopped at index ",
           i, call. = FALSE)
    }
  }

  x

}

#' @rdname MarkovChain
#' @section MarkovChainList:
#' `MarkovChainList()` constructs a container of multiple `MarkovChain` objects.
#' All elements must be of class `"MarkovChain"`. Subsetting preserves the class.
#' Subsets containing invalid elements, including `NULL` from `NA` or
#' out-of-range indices, are rejected. Empty subsets are allowed.
#'
#' @param i A numeric, character, or logical index for list subsetting. If
#'   omitted, all elements are selected.
#'
#' @return `MarkovChainList()` and its subsetting method return an object of
#'   class `"MarkovChainList"`.
#'
#' @export
MarkovChainList <- function(x) {
  validate_MarkovChainList(new_MarkovChainList(x))
}

#' @rdname MarkovChain
#' @export
`[.MarkovChainList` <- function(x, i) {
  MarkovChainList(NextMethod())
}

setOldClass(c("MarkovChainList", "list"))
