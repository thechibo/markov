# Stochastic Vector Class

Construct and validate a stochastic vector of finite, nonnegative
probabilities. The returned object has S3 class `"StochVec"`.

## Usage

``` r
new_StochVec(...)

validate_StochVec(x, tol = 1e-12)

StochVec(...)

is.StochVec(x)

# S3 method for class 'StochVec'
print(x, digits = getOption("digits"), ...)
```

## Arguments

- ...:

  For the constructors, numeric values combined into a flat vector. For
  the print method, additional arguments passed to
  [`base::print()`](https://rdrr.io/r/base/print.html).

- x:

  An object to check with `is.StochVec()`, or a `StochVec` to validate
  or print.

- tol:

  Nonnegative numeric tolerance for the absolute difference between the
  sum and one. The validator defaults to `1e-12`.

- digits:

  Number of significant digits used to round displayed values. Defaults
  to `getOption("digits")`.

## Value

The constructors and validator return a numeric vector of class
`"StochVec"`. `is.StochVec()` returns a single logical value. The print
method invisibly returns the validated vector.

## Details

If the sum differs from one by at most `1e-12`, `StochVec()` divides
every entry by that sum without a warning. A larger deviation is an
error.

## Validation

A valid `StochVec` must satisfy all of:

- all entries are finite and nonnegative

- the sum is within `1e-12` of one

- state names are unique, non-missing, and non-empty

## See also

`new_StochVec()`, `validate_StochVec()`, `is.StochVec()`

## Examples

``` r
StochVec(0.2, 0.3, 0.5)
#> Stochastic Vector
#>   length: 3 
#>     1   2   3 
#>   0.2 0.3 0.5 
x <- StochVec(c(0.1, 0.9))

# check validity / class
is.StochVec(x)
#> [1] TRUE
```
