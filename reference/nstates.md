# Number of States

Returns the number of states in a Markov chain or stochastic matrix.

## Usage

``` r
nstates(x)

# S3 method for class 'StochVec'
nstates(x)

# S3 method for class 'StochMat'
nstates(x)

# S3 method for class 'MarkovChain'
nstates(x)

# Default S3 method
nstates(x)
```

## Arguments

- x:

  A `StochVec`, a `StochMat`, or a `MarkovChain` object.

## Value

A positive integer giving the number of states.

## Examples

``` r
P <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE))

a <- StochVec(c(0.4, 0.6))
x <- MarkovChain(a = a, P = P)
nstates(a)
#> [1] 2
nstates(P)
#> [1] 2
nstates(x)
#> [1] 2
```
