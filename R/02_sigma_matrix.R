sigmaf <- function(P, n) {

  # Preliminaries
  nstate <- nrow(P)
  nr <- nstate ^ 2
  Sf <- matrix(0, nrow = nr, ncol = nr)

  # Calculate the powers of P
  list_eigen <- Pdecomp(P)
  mat_pow <- list()
  for (i in 1:(n - 1)) {
    mat_pow[[i]] <- Ppower(list_eigen, i)
  }

  for (i in 1:nstate) {
    for (j in 1:nstate) {
      for (u in 1:nstate) {
        for (v in 1:nstate) {
          k <- (i - 1) * nstate + j
          l <- (u - 1) * nstate + v
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
#' @param x
#' @param n
#'
#' @return
#' @export
#'
#' @inherit MarkovChain examples
sigma_mat <- function(x, n) {

  if (n == 1) {
    return(gamma_mat(x))
  } else {
    Sf <- sigmaf(x$P, n)
    Sigma <- Sf %*% gamma_mat(x) %*% t(Sf)
    return(Sigma)
  }
}
