# Asymptotic Covariance Matrix for k-Step Transition Probabilities

Computes the asymptotic variance-covariance matrix of the row-wise
vectorized estimator of the \\k\\-step transition matrix \\P^k\\ for an
irreducible finite-state Markov chain.

## Usage

``` r
vcov_pk(x, k, ...)

# S3 method for class 'MarkovChain'
vcov_pk(x, k, ...)

# S3 method for class 'StochMat'
vcov_pk(x, k, ...)

# Default S3 method
vcov_pk(x, k, ...)
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

- ...:

  Additional arguments passed to [`vcov_p()`](vcov_p.md).

## Value

A numeric matrix of dimension \\s^2 \times s^2\\ containing
\\\Sigma\_{\mathbf p^{(k)}}\\, where \\s\\ is the number of states.

## Details

Let \\\widehat P(n)\\ denote the empirical estimator of the one-step
transition matrix \\P\\. Then, for a fixed positive integer \\k\\, \$\$
\sqrt{n} \left\\ \widehat P^k(n) - P^k \right\\
\mathrel{\mathop{\longrightarrow}\limits^{d}} W\_{P^k}, \$\$

For which, \$\$ W\_{\mathbf p^{(k)}} =
\mathrm{vec}\left(W\_{P^k}^{\mathMk T}\right) \sim \mathcal N\_{s^2}
\left( \mathbf 0\_{s^2}, \Sigma\_{\mathbf p^{(k)}} \right). \$\$ This
function calculates the asymptotic variance covariance matrix \$\$
\Sigma\_{\mathbf p^{(k)}} = M_k \Sigma\_{\mathbf p} M_k^{\mathMk T},
\$\$

where \$\$ M_k = \sum\_{i=0}^{k-1} P^i \mathbin{\otimes}
\left(P^{k-1-i}\right)^{\mathMk T}, \$\$ and \\\otimes\\ denotes the
Kronecker product.

The one-step covariance matrix \\\Sigma\_{\mathbf p}\\ is computed by
[`vcov_p()`](vcov_p.md).

## References

GSI

A. Sadek and N. Limnios (2002). *Asymptotic properties of estimators for
nonhomogeneous Markov chains*. Statistics & Probability Letters, 60(4),
1839–1850. See Theorem 4 (p. 1847) for the asymptotic covariance matrix
of the k-step transition probabilities.

## See also

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
