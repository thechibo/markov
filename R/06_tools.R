# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tools
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Trajectory Validation            ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Validate an observed trajectory without coercing its values or attributes.
validate_trajectories <- function(x) {
  if (!is.matrix(x) && length(x) < 2) {
    stop("`x` must have length >= 2 to estimate transitions.")
  }
  if (!is.atomic(x) || (!is.null(dim(x)) && !is.matrix(x))) {
    stop("`x` must be an atomic vector or matrix of observed states.",
         call. = FALSE)
  }
  if (is.matrix(x) && (nrow(x) < 1L || ncol(x) < 2L)) {
    stop("`x` must have at least one row and two columns.", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("`x` contains NA values.", call. = FALSE)
  }
  x
}

# Common indexing for one path or matrix rows, without between-path transitions.
trajectory_data <- function(x) {
  values <- if (is.matrix(x)) as.vector(t(x)) else x
  width <- if (is.matrix(x)) ncol(x) else length(x)
  first <- seq.int(1L, length(values), by = width)
  last <- first + width - 1L
  list(values = values, first = first,
       from = seq_along(values)[-last], to = seq_along(values)[-first])
}

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
#' State names must be unique, non-missing, and non-empty.
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
  if (anyNA(names(x)) || any(names(x) == "") || anyDuplicated(names(x))) {
    stop("State names must be unique, non-missing, and non-empty.",
         call. = FALSE)
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
  if (anyNA(rownames(x)) || any(rownames(x) == "") ||
      anyDuplicated(rownames(x))) {
    stop("State names must be unique, non-missing, and non-empty.",
         call. = FALSE)
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
# Number of States                 ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Number of States
#'
#' @description
#' Returns the number of states in a Markov chain or stochastic matrix.
#'
#' @param x A `StochVec`, a `StochMat`, or a `MarkovChain` object.
#'
#' @return A positive integer giving the number of states.
#'
#' @examples
#' P <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE))
#'
#' a <- StochVec(c(0.4, 0.6))
#' x <- MarkovChain(a = a, P = P)
#' nstates(a)
#' nstates(P)
#' nstates(x)
#'
#' @name nstates
NULL

#' @rdname nstates
#' @export
nstates <- function(x) {
  UseMethod("nstates", x)
}

#' @rdname nstates
#' @export
nstates.StochVec <- function(x) {
  length(x)
}

#' @rdname nstates
#' @export
nstates.StochMat <- function(x) {
  nrow(x)
}

#' @rdname nstates
#' @export
nstates.MarkovChain <- function(x) {
  nstates(x$P)
}

#' @rdname nstates
#' @export
nstates.default <- function(x) {
  stop("`x` must be a `StochVec`, a `StochMat`, or a `MarkovChain` object.",
       call. = FALSE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Classes                          ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Communicating Classes
#'
#' @description
#' Partitions the states of a Markov chain into communicating classes.
#'
#' @param x A valid `MarkovChain` or `StochMat` object.
#' @param ... Additional arguments passed to methods. Currently unused.
#'
#' @return
#' A list of character vectors containing the state names in each
#' communicating class.
#'
#' @examples
#' P <- StochMat(matrix(c(0.5, 0.5, 0,
#'                        0.5, 0.5, 0,
#'                          0,   0, 1),
#'                      nrow = 3, byrow = TRUE))
#' mc_classes(P)
#'
#' @name mc_classes
#' @export
NULL

#' @rdname mc_classes
#' @export
mc_classes <- function(x, ...) UseMethod("mc_classes")

#' @rdname mc_classes
#' @export
mc_classes.MarkovChain <- function(x, ...) {
  x <- validate_MarkovChain(x)
  mc_classes_validated(x$P)
}

#' @rdname mc_classes
#' @export
mc_classes.StochMat <- function(x, ...) {
  mc_classes_validated(validate_StochMat(x))
}

mc_classes_validated <- function(x) {
  # Preliminaries
  state_names <- get_state_names(x)
  s <- nstates(x)
  adj <- (x > 0)
  diag(adj) <- TRUE
  id <- seq_len(s)

  reach <- adj
  for (k in seq_len(s)) {
    reach <- reach | outer(reach[id, k], reach[k, id], "&")
  }

  comm <- reach & t(reach)

  classes <- list()
  visited <- rep(FALSE, s)

  for (i in seq_len(s)) {
    if (visited[i]) {
      next
    }
    idx <- which(comm[i, ])
    visited[idx] <- TRUE
    classes[[length(classes) + 1L]] <- state_names[idx]
  }

  classes
}

#' Validate a Communicating Class
#'
#' @description
#' Checks whether a set of state names forms a communicating class.
#'
#' @param P A `StochMat` object.
#' @param c A non-empty character or integer vector containing state names.
#'
#' @return `TRUE` if `c` identifies a communicating class; otherwise `FALSE`.
#'   Invalid state names raise an error.
#'
#' @examples
#' P <- StochMat(matrix(c(0.5, 0.5, 0,
#'                        0.5, 0.5, 0,
#'                          0,   0, 1),
#'                      nrow = 3, byrow = TRUE))
#'
#' is_com_class(P, c(1L, 3L))
#' is_com_class(P, 3L)
#'
#' @export
is_com_class <- function(P, c) {
  check_state_class(P, c)$communicating
}

check_state_class <- function(P, c) {
  if (!inherits(P, "StochMat")) {
    stop("`P` must be a `StochMat` object.", call. = FALSE)
  }
  P <- validate_StochMat(P)
  if (length(c) == 0L) {
    stop("`c` must be a non-empty vector.", call. = FALSE)
  }
  if (anyNA(c)) {
    stop("`c` contains NAs.", call. = FALSE)
  }
  class_names <- as.character(c)
  state_names <- get_state_names(P)
  if (!all(class_names %in% state_names)) {
    stop("All values in `c` must be valid state names.", call. = FALSE)
  }

  classes <- mc_classes_validated(P)
  communicating <- any(vapply(classes, setequal, logical(1), y = class_names))
  list(P = P, names = class_names, communicating = communicating)
}

#' Identify a Recurrent Class
#'
#' @description
#' Checks whether a set of states forms a recurrent class.
#'
#' @param P A `StochMat` object.
#' @param c A non-empty character or integer vector containing state names.
#'
#' @return `TRUE` for a closed communicating class and `FALSE` for a transient
#'   communicating class. An error is raised if `c` does not identify a
#'   communicating class.
#'
#' @examples
#' P <- StochMat(matrix(c(0.5, 0.5, 0,
#'                        0.5, 0.5, 0,
#'                          0,   0, 1),
#'                      nrow = 3, byrow = TRUE))
#'
#' is_rec_class(P, c(1L, 2L))
#' is_rec_class(P, 3L)
#'
#' @export
is_rec_class <- function(P, c) {
  is_recurrent_class_checked(check_state_class(P, c))
}

is_recurrent_class_checked <- function(state_class) {
  if (!state_class$communicating) {
    stop("`c` must identify a communicating class.", call. = FALSE)
  }
  is_closed_class(state_class$P, state_class$names)
}

# Check a known communicating class for transitions to outside states.
is_closed_class <- function(P, class_names) {
  outside <- setdiff(get_state_names(P), class_names)
  if (length(outside) == 0L) {
    return(TRUE)
  }
  !any(P[class_names, outside, drop = FALSE] > 0)
}

#' Stationary Distribution on a Recurrent Class
#'
#' @description
#' Computes the stationary distribution associated with a recurrent
#' communicating class and embeds it in the full state space.
#'
#' @param P A `StochMat` object.
#' @param c A non-empty character or integer vector containing the state names
#'   in a recurrent class.
#'
#' @details
#' The stationary probabilities within the class are computed with the
#' Grassmann, Taksar and Heyman (1985) algorithm, using off-diagonal
#' transition probabilities.
#'
#' @return
#' A `StochVec` containing the stationary probabilities. States outside the
#' recurrent class have probability zero.
#'
#' @inherit stat references
#'
#' @importFrom stats setNames
#'
#' @examples
#' P <- StochMat(matrix(c(0.8, 0.2, 0,
#'                        0.3, 0.7, 0,
#'                          0,   0, 1),
#'                      nrow = 3, byrow = TRUE))
#'
#' stat_on_class(P, c(1L, 2L))
#' stat_on_class(P, 3L)
#'
#' @export
stat_on_class <- function(P, c) {
  state_class <- check_state_class(P, c)
  if (!is_recurrent_class_checked(state_class)) {
    stop("`c` must identify a recurrent class.", call. = FALSE)
  }
  stat_on_closed_class(state_class$P, state_class$names)
}

stat_on_closed_class <- function(P, class_names) {
  state_names <- get_state_names(P)
  Pc <- P[class_names, class_names, drop = FALSE]
  pic <- gth_stationary(Pc)
  if (max(abs(as.numeric(pic %*% Pc) - pic)) >
      sqrt(.Machine$double.eps)) {
    stop("Could not compute a valid stationary distribution.", call. = FALSE)
  }

  # Embed in full state space
  pi_full <- setNames(
    numeric(length(state_names)),
    state_names
  )
  pi_full[class_names] <- pic
  StochVec(pi_full)
}

# Grassmann, Taksar and Heyman (1985), doi:10.1287/opre.33.5.1107.
# P is the transition matrix of one closed communicating class.
gth_stationary <- function(P) {
  n <- nrow(P)
  remaining <- seq_len(n)
  eliminated <- integer(n - 1L)
  diag(P) <- 0

  for (step in seq_along(eliminated)) {
    exit_rates <- rowSums(P[remaining, remaining, drop = FALSE])
    # The largest exit rate keeps every elimination multiplier at most one.
    pivot <- which.max(exit_rates)
    state <- remaining[pivot]
    exit_rate <- exit_rates[pivot]
    if (!is.finite(exit_rate) || exit_rate <= 0) {
      stop("Could not compute a valid stationary distribution.", call. = FALSE)
    }
    eliminated[step] <- state
    remaining <- remaining[-pivot]

    incoming <- P[remaining, state] / exit_rate
    reduced <- P[remaining, remaining, drop = FALSE] +
      outer(incoming, P[state, remaining])
    diag(reduced) <- 0
    P[remaining, remaining] <- reduced
    P[remaining, state] <- incoming
  }

  # Restore the eliminated states in reverse order. Rescale at each step
  # so the weights stay bounded even when stationary probabilities vary widely.
  probabilities <- numeric(n)
  probabilities[remaining] <- 1
  for (state in rev(eliminated)) {
    weight <- sum(probabilities[remaining] * P[remaining, state])
    probabilities[remaining] <- probabilities[remaining] / (1 + weight)
    probabilities[state] <- weight / (1 + weight)
    remaining <- c(remaining, state)
  }
  probabilities
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Irreducibility                   ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Test Whether a Markov Chain is Irreducible
#'
#' @description
#' Determines whether all states of a finite Markov chain communicate.
#'
#' @param x A valid `MarkovChain` or `StochMat` object.
#' @param ... Additional arguments passed to methods. Currently unused.
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
  length(mc_classes.MarkovChain(x, ...)) == 1L
}

#' @rdname is_irreducible
#' @export
is_irreducible.StochMat <- function(x, ...) {
  length(mc_classes.StochMat(x, ...)) == 1L
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
# Matrix Powers                    ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Eigendecomposition of a Transition Matrix
#'
#' @description
#' Computes the eigendecomposition of the transpose of a transition matrix.
#'
#' @param P A square numeric transition matrix.
#'
#' @return
#' An object of class `"eigen"` containing the eigenvalues and eigenvectors
#' of `t(P)`.
#'
#' @details
#' The transpose is decomposed so that the resulting eigenvectors can be used
#' to reconstruct powers of the original transition matrix with [Ppower()].
#'
#' No validation is performed to confirm that `P` is stochastic or square.
#' Invalid inputs are handled by [base::eigen()].
#'
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
#' decomposition <- Pdecomp(P)
#' decomposition$values
#'
#' @seealso [base::eigen()], [markov::Ppower()]
#'
#' @export
Pdecomp <- function(P) {
  eigen(t(P))
}

#' Power of a Transition Matrix
#'
#' @description
#' Raises a transition matrix to a specified power using its eigendecomposition.
#'
#' @param x A numeric matrix or an `"eigen"` object.
#' @param ... Additional arguments passed to methods.
#' @param k A numeric scalar giving the exponent. For Markov-chain transition
#'   matrices, this is normally a non-negative integer.
#'
#' @return
#' A numeric or complex matrix representing the transition matrix raised to
#' the power `k`.
#'
#' @details
#' For a matrix input, `Ppower()` first computes the eigendecomposition of
#' `t(x)` using [Pdecomp()]. If
#'
#' \deqn{P^\mathsf{T} = V D V^{-1},}
#'
#' then the result is calculated as
#'
#' \deqn{P^k = \left(V D^k V^{-1}\right)^\mathsf{T}.}
#'
#' This method requires the eigenvector matrix to be invertible. An error is
#' raised by [base::solve()] when the matrix cannot be inverted. Complex
#' results may be returned when the eigendecomposition contains complex
#' eigenvalues or eigenvectors.
#'
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
#'
#' Ppower(P, k = 0)
#' Ppower(P, k = 2)
#'
#' Peigen <- Pdecomp(P)
#' Ppower(Peigen, k = 2)
#'
#' @seealso [Pdecomp()], [Ppowers()]
#'
#' @name Ppower
NULL

#' @rdname Ppower
#' @export
Ppower <- function (x, k, ...) {
  UseMethod("Ppower", x)
}

#' @rdname Ppower
#' @export
Ppower.eigen <- function(x, k, ...) {
  t(x$vectors %*% diag(x$values ^ k) %*% solve(x$vectors))
}

#' @rdname Ppower
#' @export
Ppower.matrix <- function(x, k, ...) {
  Ppower(Pdecomp(x), k = k, ...)
}

#' Successive Powers of a Transition Matrix
#'
#' @description
#' Computes successive powers of a transition matrix by repeated matrix
#' multiplication.
#'
#' @param P A nonempty square numeric matrix with finite entries.
#' @param k A single positive integer giving the highest matrix power to compute.
#'
#' @return
#' A list of `k + 1` matrices. Element `j + 1` contains \eqn{P^j}, beginning
#' with the identity matrix \eqn{P^0}. The matrices preserve the row and
#' column names of `P`.
#'
#' @details
#' Unlike [Ppower()], this function does not use an eigendecomposition.
#' Instead, each power is obtained by multiplying the preceding power by `P`:
#'
#' \deqn{P^j = P^{j-1}P.}
#'
#' Repeated multiplication can be useful when all powers from \eqn{P^0}
#' through \eqn{P^k} are required.
#'
#' The current implementation expects `k` to be at least one.
#'
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
#'
#' powers <- Ppowers(P, k = 3)
#'
#' powers[[1]] # P^0
#' powers[[2]] # P^1
#' powers[[4]] # P^3
#'
#' @seealso [Ppower()]
#'
#' @export
Ppowers <- function(P, k) {
  if (!is.matrix(P) || !is.numeric(P) || nrow(P) != ncol(P) ||
      nrow(P) == 0L || any(!is.finite(P))) {
    stop("`P` must be a nonempty square numeric matrix with finite entries.",
         call. = FALSE)
  }
  if (length(k) != 1L || !is.null(dim(k)) || !is_positive_integer(k)) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }

  n_powers <- k + 1
  Pj <- vector("list", n_powers)
  s <- nrow(P)
  Pj[[1L]] <- diag(1, s)              # P ^ 0
  dimnames(Pj[[1L]]) <- dimnames(P)
  for (j in seq.int(2L, n_powers)) {
    Pj[[j]] <- Pj[[j - 1L]] %*% P   # P ^ (j - 1)
  }
  Pj
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Other Helpers                    ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Matrix Kronecker Product
#'
#' @description
#' Computes the Kronecker product of two matrices.
#'
#' @param x A numeric matrix.
#' @param y A numeric matrix. By default, the Kronecker product of `x` with
#'   itself is computed.
#'
#' @return
#' A numeric matrix with
#' `nrow(x) * nrow(y)` rows and `ncol(x) * ncol(y)` columns.
#'
#' @details
#' The Kronecker product replaces each element of `x` with that element
#' multiplied by the complete matrix `y`.
#'
#' The function is a wrapper for [fastmatrix::kronecker.prod()].
#'
#' @examples
#' A <- matrix(c(1, 2, 3, 4), nrow = 2)
#' B <- matrix(c(0, 1, 1, 0), nrow = 2)
#'
#' kp(A, B)
#' kp(A)
#'
#' @seealso [fastmatrix::kronecker.prod()]
#'
#' @export
kp <- function(x, y = x) {
  fastmatrix::kronecker.prod(x, y)
}

#' Matrix Transformed Cross-Product
#'
#' @description
#' Computes the matrix product \eqn{C C^\mathsf{T}}, where
#' \eqn{C = A \, \mathrm{mat\_gamma}(B)}.
#'
#' @param A A numeric matrix used as the left transformation matrix.
#' @param B An object accepted by `vcov_p()`.
#'
#' @return
#' A symmetric numeric matrix with `nrow(A)` rows and columns.
#'
#' @details
#' This function is equivalent to `tcrossprod(A %*% vcov_p(B), A)`
#' provided that the matrix dimensions are conformable.
#'
#' More precisely, because the second argument to [base::tcrossprod()] is
#' supplied, the returned value is
#' \eqn{(A \, \mathrm{vcov_p}(B)) A^\mathsf{T}}.
#'
#' @examples
#' \dontrun{
#' A <- diag(2)
#' B <- c(0.5, 0.5)
#' cp(A, B)
#' }
#'
#' @seealso [base::tcrossprod()]
cp <- function(A, B) {
  tcrossprod(A %*% vcov_p(B), A)
}

#' Kronecker Delta Vector
#'
#' @description
#' Creates a vector of length `n` with a one at position `i` and zeros
#' elsewhere.
#'
#' @param i Position of the nonzero element.
#' @param n Length of the vector.
#'
#' @return A numeric vector of length `n`.
#'
#' @examples
#' cron(2, 4)
#'
#' @export
cron <- function(i, n) {
  x <- rep(0, n)
  x[i] <- 1
  x
}

#' Validate or Coerce Positive Integer Values
#'
#' @description
#' Utilities for validating and coercing positive integer-valued numeric
#' vectors.
#'
#' `is_positive_integer()` tests whether every element of `m` is finite,
#' integer-valued, greater than or equal to one, and representable as a
#' base R integer.
#'
#' `as_positive_integer()` performs the same validation and converts `m` to
#' an integer vector. It raises an error when validation fails.
#'
#' @param m An object to validate or coerce.
#'
#' @details
#' A value is considered a positive integer when it:
#'
#' - is numeric;
#' - has positive length;
#' - contains no missing or infinite values;
#' - contains only whole numbers;
#' - is greater than or equal to one; and
#' - does not exceed `.Machine$integer.max`.
#'
#' Zero, negative numbers, fractional values, `NA`, `NaN`, `Inf`, and
#' non-numeric objects are rejected.
#'
#' @return
#' `is_positive_integer()` returns a single logical value.
#'
#' `as_positive_integer()` returns `m` as an integer vector, or raises an
#' error if `m` is not a valid positive integer vector.
#'
#' @examples
#' is_positive_integer(3)
#' is_positive_integer(c(1, 2, 3))
#' is_positive_integer(c(1, 2.5))
#' is_positive_integer(0)
#'
#' as_positive_integer(3)
#' as_positive_integer(c(1, 2, 3))
#'
#' @rdname positive_integer
#' @export
is_positive_integer <- function(m) {
  !(any(!is.numeric(m)) ||
    any(!is.finite(m)) ||
    any(m < 1) ||
    any(m %% 1 != 0) ||
    any(m > .Machine$integer.max))
}

#' @rdname positive_integer
#' @export
as_positive_integer <- function(m) {
  if (!is_positive_integer(m)) {
    stop(
      sprintf(
        "`%s` must be a positive integer vector.",
        deparse1(substitute(m))
      ),
      call. = FALSE
    )
  }
  as.integer(m)
}

# Construct I - P[states, states] for a validated stochastic matrix and
# distinct character state names, preserving their supplied order.
# Sum all off-diagonal probabilities, including transitions outside the
# selected block, to avoid cancellation in 1 - P[i, i].
I_minus_P <- function(P, states) {
  off_diagonal <- P[states, , drop = FALSE]
  diagonal_indices <- cbind(seq_along(states), match(states, colnames(P)))
  off_diagonal[diagonal_indices] <- 0

  A <- -P[states, states, drop = FALSE]
  diag(A) <- rowSums(off_diagonal)
  A
}
