# Estimate a Markov Chain from One or More Trajectories

Estimates a time-homogeneous, discrete-time Markov chain from an
observed sequence of states or independent trajectories sharing a
transition matrix.

## Usage

``` r
fit(x, init = "realized")
```

## Arguments

- x:

  An atomic vector or factor of observed states in time order, or an
  atomic matrix with one independent trajectory per row and observation
  times in columns. Each trajectory must have at least two observations;
  matrices must have at least one row. Missing values are not allowed.
  State labels must be non-empty.

- init:

  Initial distribution to attach to the fitted Markov chain. Must be
  exactly one of:

  `"uniform"`

  :   Uniform distribution over the observed state space.

  `"realized"`

  :   Empirical distribution of the initial states across paths; for one
      path, concentrated on its first state. This is the default choice.

  `"stationary"`

  :   Stationary distribution implied by the estimated transition matrix
      (via [`stat()`](stat.md)).

## Value

A `MarkovChain` object.

## Details

The state space consists of the distinct observed values, reading matrix
rows in order. If `x` is numeric or ordered, the states are sorted;
otherwise, the first-appearance order is used. Numeric state labels use
additional precision when needed to distinguish values with identical
default character representations.

The function computes the empirical transition frequencies and
normalizes by row to obtain a transition probability matrix. Counts are
pooled over paths; no transitions are counted between the end of one
path and the next path.

If a state never appears as an origin (i.e., it is never observed at
times `t = 1, ..., n-1`), the corresponding transition row has no
empirical counts. In that case the function treats that state as
absorbing (probability 1 on the diagonal) and emits a warning.

The row-normalized transition matrix is validated by
[`StochMat()`](StochMat.md).

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
