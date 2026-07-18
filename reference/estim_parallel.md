# Parallel Estimation

Simulate `k` trajectories of length `m` and apply `fun`.

## Usage

``` r
estim_parallel(x, fun, k = 1000, m = 500, ...)
```

## Arguments

- x:

  an object of class `MarkovChain`.

- fun:

  the function to apply on each estimated MarkovChain.

- k:

  integer. The number of trajectories.

- m:

  integer. The length of trajectories.

- ...:

  extra arguments.

## Value

A vector holding the combined return values of `fun` for each chain.

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
pi_vec <- stat(x)

# Simulate a single trajectory of length m
m <- 1e4
traj <- sim(m, x)

# Estimate the transition probability matrix
y <- fit(traj)

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
