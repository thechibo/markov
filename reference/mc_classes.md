# Communicating Classes

Partitions the states of a Markov chain into communicating classes.

## Usage

``` r
mc_classes(x, ...)

# S3 method for class 'MarkovChain'
mc_classes(x, ...)

# S3 method for class 'StochMat'
mc_classes(x, ...)
```

## Arguments

- x:

  A `MarkovChain` or `StochMat` object.

- ...:

  Additional arguments passed to methods.

## Value

A list of character vectors containing the state names in each
communicating class.

## Examples

``` r
P <- StochMat(matrix(c(0.5, 0.5, 0,
                       0.5, 0.5, 0,
                         0,   0, 1),
                     nrow = 3, byrow = TRUE))
mc_classes(P)
#> [[1]]
#> [1] "1" "2"
#> 
#> [[2]]
#> [1] "3"
#> 
```
