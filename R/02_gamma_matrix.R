# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Gamma Matrix
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lambda_block <- function(i, P) {

  nstate <- nrow(P)
  p <- P[i, ]
  pmat <- matrix(p, nrow = nstate, ncol = 1)

  diag(p) - pmat %*% t(pmat)

}

# Sadek et Limnios (page 1846, matrix Gamma)
#' Transition Probability Vector Avar Matrix
#'
#' @param x an object of class `MarkovChain`.
#'
#' @return matrix. The asymptotic variance-covariance matrix of the transition
#' probability vector (vectorized matrix).
#' @export
#'
#' @inherit MarkovChain examples
gamma_mat <- function(x) {

  stat <- stationary(x)
  nstate <- nrow(x$P)
  list_mat <- list()

  for (i in 1:nstate){
    list_mat[[i]] <- lambda_block(i, x$P) / (stat[i] + 1e-10)
  }

  as.matrix(Matrix::bdiag(list_mat))

}
