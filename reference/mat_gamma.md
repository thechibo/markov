# Asymptotic Variance–Covariance of the Transition Probability Vector

Constructs the block-diagonal asymptotic variance–covariance matrix of
the vectorized transition probability matrix of an irreducible
finite-state Markov chain.

## Usage

``` r
mat_gamma(x, ...)

# S3 method for class 'MarkovChain'
mat_gamma(x, ...)

# S3 method for class 'StochMat'
mat_gamma(x, ...)

# Default S3 method
mat_gamma(x, ...)
```

## Arguments

- x:

  A `MarkovChain` or `StochMat` object. It must be irreducible.

- ...:

  Additional arguments passed to methods.

## Value

A numeric block-diagonal matrix of size \\s^2 \times s^2\\, where \\s\\
is the number of states.

## Details

For row \\i\\, the block \\\Gamma_i / \pi_i\\ is placed on the diagonal,
where \\\Gamma_i = \mathrm{diag}(p\_{i\bullet}) - p\_{i\bullet}
p\_{i\bullet}^\top\\ and \\\pi_i\\ is the stationary probability of
state \\i\\.

## References

A. Sadek and N. Limnios (2002). *Asymptotic properties of estimators for
nonhomogeneous Markov chains*. Statistics & Probability Letters, 60(4),
1839–1850. (See p.1846 for Γ.)

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
