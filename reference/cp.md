# Matrix Transformed Cross-Product

Computes the matrix product \\C C^\mathsf{T}\\, where \\C = A \\
\mathrm{mat\\gamma}(B)\\.

## Usage

``` r
cp(A, B)
```

## Arguments

- A:

  A numeric matrix used as the left transformation matrix.

- B:

  An object accepted by [`vcov_p()`](vcov_p.md).

## Value

A symmetric numeric matrix with `nrow(A)` rows and columns.

## Details

This function is equivalent to `tcrossprod(A %*% vcov_p(B), A)` provided
that the matrix dimensions are conformable.

More precisely, because the second argument to
[`base::tcrossprod()`](https://rdrr.io/r/base/crossprod.html) is
supplied, the returned value is \\(A \\ \mathrm{vcov_p}(B))
A^\mathsf{T}\\.

## See also

[`base::tcrossprod()`](https://rdrr.io/r/base/crossprod.html)

## Examples

``` r
if (FALSE) { # \dontrun{
A <- diag(2)
B <- c(0.5, 0.5)
cp(A, B)
} # }
```
