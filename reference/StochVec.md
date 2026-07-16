# Stochastic Vector Class

Construct and validate a stochastic vector: a numeric vector whose
components are strictly non-missing, strictly positive, and sum to
exactly 1. The returned object has S3 class `"StochVec"` and underlying
storage `"numeric"`.

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

  Numeric values coerceable to a flat numeric vector.

## Value

A `numeric` vector with additional class `"StochVec"`. If validation
fails, an error is thrown.

## Validation

A valid `StochVec` must satisfy all of:

- no `NA`s

- all entries are `> 0`

- `sum(x) == 1`

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
