power_gsi <- function(n, mylist, mat) {

  nstate <- nrow(mylist[[1]])
  s <- matrix(0, nrow = nstate ^ 2, ncol = nstate ^ 2)

  for (i in 1:n) {
    tmp <- fastmatrix::kronecker.prod(mylist[[i]], t(mylist[[n - i + 1]]))
    s <- s + tmp
  }

  s %*% mat %*% t(s)

}

