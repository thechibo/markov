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

Sadek, A., and Limnios, N. (2002). *Asymptotic properties for maximum
likelihood estimators for reliability and failure rates of Markov
chains*. Communications in Statistics - Theory and Methods, 31(10),
1837–1861.
[doi:10.1081/STA-120014916](https://doi.org/10.1081/STA-120014916) . See
Theorem 4 (p. 1847) for the component-wise covariance construction.

## See also

[`vcov_pk()`](vcov_pk.md) for the primary implementation and
[`vcov_p()`](vcov_p.md) for the one-step asymptotic covariance matrix.

## Examples

``` r
P <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE))
all.equal(markov:::vcov_pk_sl(P, k = 2L), vcov_pk(P, k = 2L))
#> [1] TRUE
```
