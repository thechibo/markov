# Stationary Distribution of a Markov Chain

Computes stationary distribution(s) of a finite-state, time-homogeneous
Markov chain.

## Usage

``` r
stat(x, ...)
```

## Arguments

- x:

  A valid `MarkovChain` or `StochMat` object.

- ...:

  Additional arguments passed to methods. Currently unused.

## Value

A `StochVec` (unique recurrent class), or a list of `StochVec` objects.

## Details

For a transition matrix \\P\\, a stationary distribution \\\pi\\
satisfies \$\$\pi^\top P = \pi^\top,\qquad \sum_i \pi_i = 1,\\ \pi_i \ge
0.\$\$ For a `MarkovChain`, both components are validated and the
stationary distribution depends only on the transition matrix.

If the chain has a single closed communicating class, the stationary
distribution is unique (though the chain may be periodic). If there are
multiple closed communicating classes, stationary distributions are not
unique. Every stationary distribution is a convex combination of the
extremal ones, each supported on one closed class.

The Grassmann, Taksar and Heyman (1985) algorithm, also known as GTH,
computes the stationary probabilities within each closed class using
off-diagonal transitions.

This function returns:

- a single `StochVec` if the stationary distribution is unique;

- otherwise a list of `StochVec` objects, containing extremal stationary
  distributions of the corresponding recurrent classes.

## References

Grassmann, W. K., Taksar, M. I., and Heyman, D. P. (1985). Regenerative
Analysis and Steady State Distributions for Markov Chains. *Operations
Research*, 33(5), 1107–1116.
[doi:10.1287/opre.33.5.1107](https://doi.org/10.1287/opre.33.5.1107)

## See also

[`MarkovChain()`](MarkovChain.md), [`sim()`](sim.md), [`fit()`](fit.md)

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
