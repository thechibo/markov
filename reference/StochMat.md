# Stochastic Matrix Class

Construct and validate a stochastic matrix: a square numeric matrix
whose entries are strictly non-missing, strictly positive, and each row
sums to exactly 1. The returned object has S3 class `"StochMat"` and
underlying storage `"matrix"`.

## Usage

``` r
new_StochMat(x, ...)

validate_StochMat(x, tol = 1e-12)

StochMat(...)

is.StochMat(x)

# S3 method for class 'StochMat'
print(x, digits = getOption("digits"), ...)
```

## Arguments

- x:

  For `new_StochMat()`, a matrix or something coercible to one.

- ...:

  Objects coercible to a matrix. Passed to
  [`matrix()`](https://rdrr.io/r/base/matrix.html) when needed.

## Value

A matrix with additional class `"StochMat"`. If validation fails, an
error is thrown.

## Validation

A valid `StochMat` must satisfy all of:

- is a square matrix (`nrow(x) == ncol(x)`)

- no `NA`s

- all entries are `> 0`

- `rowSums(x) == 1` for every row

## See also

`new_StochMat()`, `validate_StochMat()`, `is.StochMat()`

## Examples

``` r
StochMat(diag(3))
#> Stochastic Matrix
#>   dim: 3 x 3 
#>     1 2 3
#>   1 1 0 0
#>   2 0 1 0
#>   3 0 0 1
StochMat(matrix(c(0.2,0.8, 0.7,0.3), ncol = 2, byrow = TRUE))
#> Stochastic Matrix
#>   dim: 2 x 2 
#>       1   2
#>   1 0.2 0.8
#>   2 0.7 0.3

m <- StochMat(diag(2))
is.StochMat(m)
#> [1] TRUE
```
