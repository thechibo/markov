# Simulate a Markov Chain Trajectory

Generate a single sample path of length `m` from a finite-state,
time-homogeneous Markov chain.

## Usage

``` r
mc_sim(m, x)
```

## Arguments

- m:

  Positive integer. The desired trajectory length.

- x:

  A `MarkovChain` object.

## Value

A factor of length `m` giving the visited states in order. Levels
correspond to the state labels (from `names(x$a)` or `rownames(x$P)` if
available).

## Details

The initial state is drawn from the chain's initial distribution `a`,
and subsequent states are generated iteratively using the transition
matrix `P`.

## See also

[`MarkovChain()`](MarkovChain.md), [`mc_estim()`](mc_estim.md),
[`mc_stat()`](mc_stat.md)

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
