# Identify a Recurrent Class

Checks whether a set of states forms a recurrent class.

## Usage

``` r
is_rec_class(P, c)
```

## Arguments

- P:

  A `StochMat` object.

- c:

  A non-empty character or integer vector containing state names.

## Value

`TRUE` for a closed communicating class and `FALSE` for a transient
communicating class. An error is raised if `c` does not identify a
communicating class.

## Examples

``` r
P <- StochMat(matrix(c(0.5, 0.5, 0,
                       0.5, 0.5, 0,
                         0,   0, 1),
                     nrow = 3, byrow = TRUE))

is_rec_class(P, c(1L, 2L))
#> [1] TRUE
is_rec_class(P, 3L)
#> [1] TRUE
```
