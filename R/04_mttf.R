# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mean Time to Failure                                                      ----
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
mttf <- function(x, up) {

  r <- length(up)
  vec_ones <- rep(1, r)
  B <- solve(diag(r) - x$P[up, up])

  as.numeric(x$a[up] %*% B %*% vec_ones)

}

#' @rdname mttf
#' @export
var_mttf <- function(x, up) {

  s <- nstates(x)
  r <- length(up)
  vec_ones <- rep(1, r)
  B <- solve(diag(r) - x$P[up, up])

  M <- kp(x$a[up] %*% B, t(vec_ones) %*% t(B))

  upmat <- matrix(FALSE, nrow = s, ncol = s)
  upmat[up, up] <- TRUE
  upvec <- as.vector(upmat)

  as.vector(M %*% mat_gamma(x)[upvec, upvec] %*% t(M))

}

# Sadek and Limnios, Theorem 8
#' @rdname mttf
#' @export
var_mttf_sl <- function(x, up) {
  psi_prime <- deriv_psi(x, up)

  upmat <- matrix(FALSE, nrow = nstates(x), ncol = nstates(x))
  upmat[up, up] <- TRUE
  upvec <- as.vector(upmat)

  as.vector(t(psi_prime) %*% mat_gamma(x)[upvec, upvec] %*% psi_prime)

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
