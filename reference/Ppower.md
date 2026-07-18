# Power of a Transition Matrix

Raises a transition matrix to a specified power using its
eigendecomposition.

## Usage

``` r
Ppower(x, k, ...)

# S3 method for class 'eigen'
Ppower(x, k, ...)

# S3 method for class 'matrix'
Ppower(x, k, ...)
```

## Arguments

- x:

  A numeric matrix or an `"eigen"` object.

- k:

  A numeric scalar giving the exponent. For Markov-chain transition
  matrices, this is normally a non-negative integer.

- ...:

  Additional arguments passed to methods.

## Value

A numeric or complex matrix representing the transition matrix raised to
the power `k`.

## Details

For a matrix input, `Ppower()` first computes the eigendecomposition of
`t(x)` using [`Pdecomp()`](Pdecomp.md). If

\$\$P^\mathsf{T} = V D V^{-1},\$\$

then the result is calculated as

\$\$P^k = \left(V D^k V^{-1}\right)^\mathsf{T}.\$\$

This method requires the eigenvector matrix to be invertible. An error
is raised by [`base::solve()`](https://rdrr.io/r/base/solve.html) when
the matrix cannot be inverted. Complex results may be returned when the
eigendecomposition contains complex eigenvalues or eigenvectors.

## See also

[`Pdecomp()`](Pdecomp.md), [`Ppowers()`](Ppowers.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)

Ppower(P, k = 0)
#>              [,1]          [,2]
#> [1,] 1.000000e+00 -1.350959e-16
#> [2,] 3.707416e-17  1.000000e+00
Ppower(P, k = 2)
#>      [,1] [,2]
#> [1,] 0.70 0.30
#> [2,] 0.45 0.55

Peigen <- Pdecomp(P)
Ppower(Peigen, k = 2)
#>      [,1] [,2]
#> [1,] 0.70 0.30
#> [2,] 0.45 0.55
```
