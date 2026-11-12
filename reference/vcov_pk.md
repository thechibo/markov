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

  Additional arguments passed to methods. Currently unused.

## Value

A numeric matrix of dimension \\s^2 \times s^2\\ containing
\\\Sigma\_{\mathbf p^{(k)}}\\, where \\s\\ is the number of states.

## Details

Let \\\widehat P(n)\\ denote the empirical estimator of the one-step
transition matrix \\P\\. Then, for a fixed positive integer \\k\\, \$\$
\sqrt{n} \left\\ \widehat P^k(n) - P^k \right\\
\mathrel{\mathop{\longrightarrow}\limits^{d}} W\_{P^k}, \$\$

For which, \$\$ W\_{\mathbf p^{(k)}} =
\mathrm{vec}\left(W\_{P^k}^{\top}\right) \sim \mathcal N\_{s^2} \left(
\mathbf 0\_{s^2}, \Sigma\_{\mathbf p^{(k)}} \right). \$\$ This function
calculates the asymptotic variance covariance matrix \$\$
\Sigma\_{\mathbf p^{(k)}} = M_k \Sigma\_{\mathbf p} M_k^{\top}, \$\$

where \$\$ M_k = \sum\_{i=0}^{k-1} P^i \mathbin{\otimes}
\left(P^{k-1-i}\right)^{\top}, \$\$ and \\\otimes\\ denotes the
Kronecker product.

The one-step covariance matrix \\\Sigma\_{\mathbf p}\\ is computed by
[`vcov_p()`](vcov_p.md).

## References

Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026). *Matrix
asymptotic calculus for plug-in maximum likelihood estimators in finite
Markov chains*. arXiv preprint arXiv:2607.17161.
<https://arxiv.org/abs/2607.17161>. See Proposition 4 for the matrix
differential and its Kronecker representation.

Sadek, A., and Limnios, N. (2002). *Asymptotic properties for maximum
likelihood estimators for reliability and failure rates of Markov
chains*. Communications in Statistics - Theory and Methods, 31(10),
1837–1861.
[doi:10.1081/STA-120014916](https://doi.org/10.1081/STA-120014916) . See
Theorem 4 (p. 1847) for the asymptotic covariance matrix of the k-step
transition probabilities.

## See also

[`vcov_p()`](vcov_p.md) for the one-step asymptotic covariance matrix.

## Examples

``` r
P <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE))
dim(vcov_pk(P, k = 2L))
#> [1] 4 4
identical(vcov_pk(P, k = 1L), vcov_p(P))
#> [1] TRUE
```
