# Markov Chain (S3) Class

Construct and validate a finite‐state, discrete‐time Markov chain,
represented by an initial distribution vector `a` and a row‐stochastic
transition probability matrix `P`. Both are validated and coerced to
classes [`StochVec`](StochVec.md) and [`StochMat`](StochMat.md).

## Usage

``` r
new_MarkovChain(a, P)

validate_MarkovChain(x)

MarkovChain(a, P)

is.MarkovChain(x)

# S3 method for class 'MarkovChain'
print(x, digits = getOption("digits"), ...)

new_MarkovChainList(x)

validate_MarkovChainList(x)

MarkovChainList(x)

# S3 method for class 'MarkovChainList'
x[i]
```

## Arguments

- a:

  Numeric, initial law.

- P:

  Numeric matrix, transition probabilities (square row‐stochastic).

- x:

  For `MarkovChainList()`, a list of `MarkovChain` objects.

## Value

A list with class `"MarkovChain"` and elements `a` and `P`. Errors on
validation failure.

For `MarkovChainList()`, an object of class `"MarkovChainList"`.

## Structure

A `MarkovChain` is a list with two named elements:

- a:

  initial distribution, of class `StochVec`

- P:

  transition matrix, of class `StochMat`

## Validation

A valid `MarkovChain` must satisfy:

- has exactly two elements named `"a"` and `"P"`

- `a` validates as `StochVec`

- `P` validates as `StochMat`

## MarkovChainList

`MarkovChainList()` constructs a container of multiple `MarkovChain`
objects. All elements must be of class `"MarkovChain"`. Subsetting
preserves the class.

## See also

[`StochVec()`](StochVec.md), [`StochMat()`](StochMat.md),
`is.MarkovChain()`

## Examples

``` r
# Choose the transition probability matrix and the initial law of the MC
a <- c(0.5, 0.5, 0)
P <- matrix(c(0.3, 0.7, 0,
              0.6, 0.1, 0.3,
              0.5, 0,   0.5), nrow = 3, ncol = 3, byrow = TRUE)

# Create a Markov Chain
x <- MarkovChain(a = a, P = P)

# Find Stationary distribution
pi_vec <- mc_stat(x)

# Simulate a single trajectory of length m
m <- 1e4
traj <- mc_sim(m, x)

# Estimate the transition probability matrix
y <- mc_estim(traj)

# Make a list of Markov Chains
L <- MarkovChainList(list(x, y))
L[[1]]
#> Markov Chain
#>   Number of states: 3 
#> 
#>   Initial distribution (a):
#>       1   2   3 
#>     0.5 0.5 0.0 
#> 
#>   Transition matrix (P):
#>         1   2   3
#>     1 0.3 0.7 0.0
#>     2 0.6 0.1 0.3
#>     3 0.5 0.0 0.5
```
