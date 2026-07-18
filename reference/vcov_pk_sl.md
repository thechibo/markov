# Reference Implementation of the k-Step Asymptotic Covariance Matrix

Computes the asymptotic variance-covariance matrix of the row-wise
vectorized estimator of the \\k\\-step transition matrix \\P^k\\ for an
irreducible finite-state Markov chain.

## Usage

``` r
vcov_pk_sl(x, k)
```

## Arguments

- x:

  A `MarkovChain` or `StochMat` object representing an \\s \times s\\
  transition matrix \\P\\. The associated Markov chain must be
  irreducible.

- k:

  A positive integer specifying the transition step. Thus, \\k = 1\\
  corresponds to the one-step transition matrix, \\k = 2\\ to the
  two-step transition matrix, and so forth.

## Value

A dense numeric matrix of dimension \\s^2 \times s^2\\ containing the
asymptotic variance-covariance matrix \\\Sigma\_{\mathbf p(k)}\\, where
\\s\\ is the number of states.

## Details

Unlike [`vcov_pk()`](vcov_pk.md), which uses a matrix and
Kronecker-product implementation, this function follows the
corresponding scalar or component-wise expressions. It is used for
internal testing only.

## References

A. Sadek and N. Limnios (2002). *Asymptotic properties of estimators for
nonhomogeneous Markov chains*. Statistics & Probability Letters, 60(4),
1839–1850.

## See also

[`vcov_pk()`](vcov_pk.md) for the primary implementation and
[`vcov_p()`](vcov_p.md) for the one-step asymptotic covariance matrix.

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
