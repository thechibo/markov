# Test Whether a Markov Chain is Irreducible

Determines whether all states of a finite Markov chain communicate.

## Usage

``` r
is_irreducible(x, ...)

# S3 method for class 'MarkovChain'
is_irreducible(x, ...)

# S3 method for class 'StochMat'
is_irreducible(x, ...)

# Default S3 method
is_irreducible(x, ...)
```

## Arguments

- x:

  A `MarkovChain` or `StochMat` object.

- ...:

  Additional arguments passed to methods.

## Value

A single logical value: `TRUE` if the chain is irreducible and `FALSE`
otherwise.
