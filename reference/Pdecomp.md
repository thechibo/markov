# Eigendecomposition of a Transition Matrix

Computes the eigendecomposition of the transpose of a transition matrix.

## Usage

``` r
Pdecomp(P)
```

## Arguments

- P:

  A square numeric transition matrix.

## Value

An object of class `"eigen"` containing the eigenvalues and eigenvectors
of `t(P)`.

## Details

The transpose is decomposed so that the resulting eigenvectors can be
used to reconstruct powers of the original transition matrix with
[`Ppower()`](Ppower.md).

No validation is performed to confirm that `P` is stochastic or square.
Invalid inputs are handled by
[`base::eigen()`](https://rdrr.io/r/base/eigen.html).

## See also

[`base::eigen()`](https://rdrr.io/r/base/eigen.html),
[`Ppower()`](Ppower.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
decomposition <- Pdecomp(P)
decomposition$values
#> [1] 1.0 0.5
```
