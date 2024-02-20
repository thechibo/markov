# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tools
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Kronecker
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cron <- function(i, n) {
  x <- rep(0, n)
  x[i] <- 1
  x
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Matrix Powers
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pdecomp <- function(P) {
  eigen(t(P))
}

Ppower <- function (x, ...) {
  UseMethod("Ppower", x)
}

Ppower.eigen <- function(r, k) {
  rvec <- r$vectors[, which.min(abs(r$values - 1))]
  rvec / sum(rvec)
  t(r$vectors %*% diag(r$values ^ k) %*% solve(r$vectors))
}

Ppower.matrix <- function(P, k) {
  r <- Pdecomp(P)
  rvec <- r$vectors[, which.min(abs(r$values - 1))]
  rvec / sum(rvec)
  t(r$vectors %*% diag(r$values ^ k) %*% solve(r$vectors))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parallel
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`%dopar%` <- foreach::`%dopar%`
`%do%` <- foreach::`%do%`
