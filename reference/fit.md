# Estimate a Markov chain from an observed state sequence

Estimates a time-homogeneous, discrete-time Markov chain from an
observed sequence of states.

## Usage

``` r
fit(x, init = c("realized", "uniform", "stationary"))
```

## Arguments

- x:

  A vector of observed states in time order. Can be numeric, character,
  factor, or ordered. Must have length at least 2 and contain no `NA`s.

- init:

  Initial distribution to attach to the fitted Markov chain. Must be
  exactly one of:

  `"uniform"`

  :   Uniform distribution over the observed state space.

  `"realized"`

  :   Degenerate distribution concentrated on the first observed state
      `x[1]`. This is the default choice.

  `"stationary"`

  :   Stationary distribution implied by the estimated transition matrix
      (via [`stat()`](stat.md)).

## Value

A `MarkovChain` object.

## Details

The state space is defined as `unique(x)`. If `x` is numeric or ordered,
the states are sorted; otherwise, the first-appearance order is used.

The function computes the empirical transition frequencies and
normalizes by row to obtain a transition probability matrix.

If a state never appears as an origin (i.e., it is never observed at
times `t = 1, ..., n-1`), the corresponding transition row has no
empirical counts. In that case the function treats that state as
absorbing (probability 1 on the diagonal) and emits a warning.

A numerical cleanup step will rescale rows if their sums deviate from 1
by by more than 1e-12 (tolerance), emitting a warning.

## See also

[`MarkovChain()`](MarkovChain.md), [`sim()`](sim.md),
[`stat()`](stat.md)

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
