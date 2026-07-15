# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sigma Matrix                                                              ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Need to improve documentation.
## Need to improve sigmaf efficiency.

#' Sigma Matrix
#'
#' @param x an object of class `MarkovChain`.
#' @param k integer. The power.
#'
#' @return matrix.
#' @export
#'
#' @importFrom fastmatrix kronecker.prod
#'
#' @inherit MarkovChain examples
mat_sigma <- function(x, k) {

  if (k == 1) {
    return(mat_gamma(x))
  } else {
    s <- nstates(x)

    # Compute the matrix powers
    Pk <- Ppowers(P, k)
    Pkt <- lapply(Pk, t)

    Sf <- matrix(0, s ^ 2, s ^ 2)
    for (i in 1:k) {
      Sf <- Sf + fastmatrix::kronecker.prod(Pk[[i]], Pkt[[k - i + 1]])
    }

    Sf %*% mat_gamma(x) %*% t(Sf)
  }
}

sigmaf <- function(P, n) {

  # Preliminaries
  s <- nrow(P)
  nr <- s ^ 2
  Sf <- matrix(0, nrow = nr, ncol = nr)

  # Calculate the powers of P
  list_eigen <- Pdecomp(P)
  mat_pow <- list()
  for (i in 1:(n - 1)) {
    mat_pow[[i]] <- Ppower(list_eigen, i)
  }

  for (i in 1:s) {
    for (j in 1:s) {
      for (u in 1:s) {
        for (v in 1:s) {
          k <- (i - 1) * s + j
          l <- (u - 1) * s + v
          Sf[k, l] <- 1 * (j == v) * mat_pow[[n - 1]][i, u] +
                      1 * (i == u) * mat_pow[[n - 1]][v, j]
          if (n > 2) {
            s <- 0
            for (m in 2:(n - 1)) {
              s <- s + mat_pow[[n - m]][i, u] * mat_pow[[m - 1]][v, j]
            }
            Sf[k, l] <- Sf[k, l] + s
          }
        }
      }
    }
  }

  Sf

}

#' Title
#'
#' @param x an object of class `MarkovChain`.
#' @param k integer.
#'
#' @return a matrix.
#' @export
#'
#' @inherit MarkovChain examples
mat_sigma_sl <- function(x, k) {

  if (k == 1) {
    return(mat_gamma(x))
  } else {
    Sf <- sigmaf(x$P, k)
    return(Sf %*% mat_gamma(x) %*% t(Sf))
  }
}
