# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mean Time to Failure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Mean Time to Failure
#'
#' @param x an object of class `MarkovChain`.
#' @param up integer or character. The up states.
#' @param ... extra arguments.
#'
#' @return numeric. The Mean Time to Failure.
#' @export
#'
#' @inherit MarkovChain examples
mttf <- function (x, ...) {
  UseMethod("mttf", x)
}

#' @export
mttf.MarkovChain <- function(x, up) {

  # Vector of the conditional mean times to failure
  vec_ones <- rep(1, length(up))
  mat <- diag(length(up)) - x$P[up, up]

  # Mean Time to Failure
  as.numeric(x$a[up] %*% solve(mat) %*% vec_ones)

}

# Sadek and Limnios, Theorem 8
#' @rdname mttf
#' @export
var_mttf <- function(x, up) {
  psi_prime <- deriv_psi(x, up)
  t(psi_prime) %*% gamma_mat(x)[up, up] %*% psi_prime
}

#' @rdname mttf
#' @export
deriv_psi <- function(x, up) {

  s1 <- length(up)
  psi_prime <- numeric(s1 ^ 2)
  col_ones <- matrix(1, ncol = 1, nrow = s1)

  # inverse of (I - P11)
  Lambda <- solve(diag(s1) - x$P[up, up])

  for (i in 1:(s1 ^ 2)) {
    mat <- matrix(cron(i, s1 ^ 2), nrow = s1, ncol = s1)
    psi_prime[i] <- x$a[up] %*% Lambda %*% mat %*% Lambda %*% col_ones
  }

  psi_prime

}
