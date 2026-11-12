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

The implementation applies the matrix delta method using a Kronecker
product: \$\$\Sigma\_\pi = M\\\Gamma\\M^\top,\$\$ where \\\Gamma\\ is
the asymptotic covariance matrix for the row-wise vectorized transition
estimator \\\mathrm{vec}(\hat P^\top)\\ returned by
[`vcov_p()`](vcov_p.md). With \\\pi\\ written as a row vector, \$\$M =
\pi \otimes (I - P^\top + \mathbf{1}\mathbf{1}^\top)^{-1}.\$\$

## References

Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026). *Matrix
asymptotic calculus for plug-in maximum likelihood estimators in finite
Markov chains*. arXiv preprint arXiv:2607.17161.
<https://arxiv.org/abs/2607.17161>. See Proposition 7 for the stationary
differential and covariance.

## See also

[`stat()`](stat.md), [`vcov_p()`](vcov_p.md),
[`is_irreducible()`](is_irreducible.md)

## Examples

``` r
P <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE))
Sigma <- vcov_stat(P)
dim(Sigma) # One stationary probability per state
#> [1] 2 2
```
