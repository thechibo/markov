# Get or Set State Names

Utilities for assigning and retrieving the state names associated with
stochastic vectors, stochastic matrices, and Markov chains.

`set_state_names()` assigns default sequential state names when names
are missing. Existing names are preserved.

`get_state_names()` retrieves the state names attached to an object.

## Usage

``` r
set_state_names(x, ...)

# S3 method for class 'StochVec'
set_state_names(x, ...)

# S3 method for class 'StochMat'
set_state_names(x, ...)

# Default S3 method
set_state_names(x, ...)

get_state_names(x, ...)

# S3 method for class 'StochVec'
get_state_names(x, ...)

# S3 method for class 'StochMat'
get_state_names(x, ...)

# S3 method for class 'MarkovChain'
get_state_names(x, ...)

# Default S3 method
get_state_names(x, ...)
```

## Arguments

- x:

  A `StochVec`, `StochMat`, or, for `get_state_names()`, a `MarkovChain`
  object.

- ...:

  Additional arguments passed to methods. Currently unused.

## Value

`set_state_names()` returns `x` with state names assigned where
necessary.

`get_state_names()` returns a character vector containing the state
names.

## Details

For an unnamed `StochVec`, `set_state_names()` assigns the names `"1"`,
`"2"`, ..., according to the vector length.

For a `StochMat`, missing row or column names are assigned as `"1"`,
`"2"`, ..., according to the matrix dimensions. After assignment, the
row and column names must be identical. An error is raised when they
differ.

State names must be unique, non-missing, and non-empty.

For a `MarkovChain`, `get_state_names()` retrieves the row names of its
transition matrix.

## Examples

``` r
a <- StochVec(c(0.4, 0.6))
get_state_names(a)
#> [1] "1" "2"

P <- StochMat(matrix(c(0.8, 0.2,
                       0.3, 0.7),
                     nrow = 2, byrow = TRUE))
get_state_names(P)
#> [1] "1" "2"

named_a <- StochVec(c(low = 0.4, high = 0.6))
get_state_names(named_a)
#> [1] "low"  "high"

named_P <- StochMat(
  matrix(c(0.8, 0.2,
           0.3, 0.7),
         nrow = 2, byrow = TRUE,
         dimnames = list(c("low", "high"), c("low", "high")))
)
get_state_names(named_P)
#> [1] "low"  "high"

mc <- MarkovChain(a = named_a, P = named_P)
get_state_names(mc)
#> [1] "low"  "high"
```
