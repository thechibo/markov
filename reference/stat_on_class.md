# Stationary Distribution on a Recurrent Class

Computes the stationary distribution associated with a recurrent
communicating class and embeds it in the full state space.

## Usage

``` r
stat_on_class(P, c)
```

## Arguments

- P:

  A `StochMat` object.

- c:

  A non-empty character or integer vector containing the state names in
  a recurrent class.

## Value

A `StochVec` containing the stationary probabilities. States outside the
recurrent class have probability zero.

## Examples

``` r
P <- StochMat(matrix(c(0.8, 0.2, 0,
                       0.3, 0.7, 0,
                         0,   0, 1),
                     nrow = 3, byrow = TRUE))

stat_on_class(P, c(1L, 2L))
#> Stochastic Vector
#>   length: 3 
#>     1   2   3 
#>   0.6 0.4 0.0 
stat_on_class(P, 3L)
#> Stochastic Vector
#>   length: 3 
#>   1 2 3 
#>   0 0 1 
```
