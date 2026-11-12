# Stochastic Matrix Class

Construct and validate a square matrix of finite, nonnegative transition
probabilities. The returned object has S3 class `"StochMat"`.

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

  For `new_StochMat()`, a matrix or data for
  [`base::matrix()`](https://rdrr.io/r/base/matrix.html). For
  `is.StochMat()`, the object to check. For validation or printing, a
  `StochMat`.

- ...:

  For the constructors, data and arguments for
  [`base::matrix()`](https://rdrr.io/r/base/matrix.html) when the input
  is supplied as data rather than as a matrix. For the print method,
  additional arguments passed to
  [`base::print()`](https://rdrr.io/r/base/print.html).

- tol:

  Nonnegative numeric tolerance for the absolute difference between each
  row sum and one after clipping. The validator defaults to `1e-12`.

- digits:

  Number of significant digits used to round displayed values. Defaults
  to `getOption("digits")`.

## Value

The constructors and validator return a matrix of class `"StochMat"`.
`is.StochMat()` returns a single logical value. The print method
invisibly returns the validated matrix.

## Details

Input entries below `-1e-16` are rejected. Entries in `[-1e-16, 0)` are
set to zero. After clipping, each row sum must differ from one by at
most `1e-12`. A larger deviation is an error. Each accepted row is
divided by its sum. `StochMat()` warns only if it clips negative
entries.

## Validation

A valid `StochMat` must satisfy all of:

- is a square matrix (`nrow(x) == ncol(x)`)

- has at least one state

- all entries are finite and nonnegative

- every row sum is within `1e-12` of one

- state names are unique, non-missing, and non-empty

Assign the returned object to retain corrections. `is.StochMat()` only
checks validity.

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
