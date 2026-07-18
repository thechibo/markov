# Validate a Communicating Class

Checks whether a set of state names forms a communicating class.

## Usage

``` r
is_com_class(P, c)
```

## Arguments

- P:

  A `StochMat` object.

- c:

  A non-empty character or integer vector containing state names.

## Value

`TRUE` if `c` identifies a communicating class; otherwise `FALSE`.

## Examples

``` r
P <- StochMat(matrix(c(0.5, 0.5, 0,
                       0.5, 0.5, 0,
                         0,   0, 1),
                     nrow = 3, byrow = TRUE))

is_com_class(P, c(1L, 3L))
#> [1] FALSE
is_com_class(P, 3L)
#> [1] TRUE
```
