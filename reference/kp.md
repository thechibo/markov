# Matrix Kronecker Product

Computes the Kronecker product of two matrices.

## Usage

``` r
kp(x, y = x)
```

## Arguments

- x:

  A numeric matrix.

- y:

  A numeric matrix. By default, the Kronecker product of `x` with itself
  is computed.

## Value

A numeric matrix with `nrow(x) * nrow(y)` rows and `ncol(x) * ncol(y)`
columns.

## Details

The Kronecker product replaces each element of `x` with that element
multiplied by the complete matrix `y`.

The function is a wrapper for
[`fastmatrix::kronecker.prod()`](https://rdrr.io/pkg/fastmatrix/man/kronecker.prod.html).

## See also

[`fastmatrix::kronecker.prod()`](https://rdrr.io/pkg/fastmatrix/man/kronecker.prod.html)

## Examples

``` r
A <- matrix(c(1, 2, 3, 4), nrow = 2)
B <- matrix(c(0, 1, 1, 0), nrow = 2)

kp(A, B)
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    1    0    3
#> [2,]    1    0    3    0
#> [3,]    0    2    0    4
#> [4,]    2    0    4    0
kp(A)
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    3    3    9
#> [2,]    2    4    6   12
#> [3,]    2    6    4   12
#> [4,]    4    8    8   16
```
