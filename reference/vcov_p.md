# Asymptotic Variance–Covariance of the Transition Probability Vector

Constructs the block-diagonal asymptotic variance–covariance matrix for
the row-wise vectorization of an empirical transition matrix from an
irreducible finite-state Markov chain.

## Usage

``` r
vcov_p(x, ...)

# S3 method for class 'MarkovChain'
vcov_p(x, ...)

# S3 method for class 'StochMat'
vcov_p(x, ...)

# Default S3 method
vcov_p(x, ...)
```

## Arguments

- x:

  A `MarkovChain` or `StochMat` object representing an \\s \times s\\
  row-stochastic transition matrix. The Markov chain must be
  irreducible.

- ...:

  Additional arguments passed to methods.

## Value

A numeric block-diagonal matrix of dimension \\s^2 \times s^2\\, where
\\s\\ is the number of states.

## Details

Let \\\widehat{\mathbf P}(n)\\ be the empirical transition matrix
estimated from a trajectory of length \\n\\, and let \\\widehat{\mathbf
p}(n)\\ be its row-wise vectorization. Under irreducibility, the
empirical transition matrix is strongly consistent and asymptotically
normal: \$\$ \sqrt{n} \left\\ \widehat{\mathbf p}(n) - \mathbf p
\right\\ \xrightarrow\[n\to\infty\]{d} \mathbf W\_{\mathbf p}, \qquad
\mathbf W\_{\mathbf p} \sim \mathcal N\_{s^2} \left( \mathbf 0\_{s^2},
\Sigma\_{\mathbf p} \right). \$\$

This function calculates the matrix \\\Sigma\_{\mathbf p}\\.
Specifically, if \\\boldsymbol{\pi} = (\pi_1,\ldots,\pi_s)^{\mathsf T}\\
is the unique stationary distribution, then the asymptotic covariance
matrix is \$\$ \Sigma\_{\mathbf p} = \begin{pmatrix} \pi_1^{-1}\Lambda_1
& 0 & \cdots & 0 \\ 0 & \pi_2^{-1}\Lambda_2 & \cdots & 0 \\ \vdots &
\vdots & \ddots & \vdots \\ 0 & 0 & \cdots & \pi_s^{-1}\Lambda_s
\end{pmatrix}, \$\$

where \$\$ \Lambda_i = \operatorname{diag}(p\_{i\bullet}) -
p\_{i\bullet}p\_{i\bullet}^{\mathsf T}, \qquad i = 1,\ldots,s. \$\$

## References

GSI - Theorem 1

A. Sadek and N. Limnios (2002). *Asymptotic properties of estimators for
nonhomogeneous Markov chains*. Statistics & Probability Letters, 60(4),
1839–1850. See Proposition 1 (p. 1846) for the asymptotic covariance
matrix of the transition probabilities.

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
