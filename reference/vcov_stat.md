# Asymptotic Variance–Covariance Matrix of the stationary distribution MLE

Computes the asymptotic variance–covariance matrix of the maximum
likelihood estimator (MLE) of the stationary distribution of an
irreducible, finite-state Markov chain.

## Usage

``` r
vcov_stat(x, ...)

# S3 method for class 'MarkovChain'
vcov_stat(x, ...)

# S3 method for class 'StochMat'
vcov_stat(x, ...)

# Default S3 method
vcov_stat(x, ...)
```

## Arguments

- x:

  A `MarkovChain` or `StochMat` object. Must be irreducible.

- ...:

  Additional arguments passed to methods.

## Value

A numeric \\s \times s\\ matrix, where \\s\\ is the number of states,
giving the asymptotic covariance matrix \\\Sigma\_\pi\\.

## Details

Let \\P\\ be the transition matrix of an irreducible finite-state Markov
chain with stationary distribution \\\pi\\. If \\\hat P\\ is the MLE of
\\P\\ based on a long trajectory, then the induced estimator \\\hat\pi =
\pi(\hat P)\\ satisfies \$\$\sqrt{n}(\hat\pi - \pi) \xrightarrow{d}
\mathcal{N}(0, \Sigma\_\pi).\$\$

The implementation uses the delta method in matrix form and Kronecker:
\$\$\Sigma\_\pi = M\\\Gamma\\M^\top,\$\$ where \\\Gamma\\ is the
asymptotic covariance matrix of \\\mathrm{vec}(\hat P)\\ returned by
[`vcov_p()`](vcov_p.md), and \$\$M = \pi^\top \otimes (I - P^\top +
\mathbf{1}\mathbf{1}^\top)^{-1}.\$\$

## References

GSI

## See also

[`stat()`](stat.md), [`vcov_p()`](vcov_p.md),
[`is_irreducible()`](is_irreducible.md)

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
