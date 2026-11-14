# Component-wise Reference for Reliability Variance

Compute the pointwise asymptotic variance using the component-wise
covariance expansion of Theorem 6 of Sadek and Limnios (2002, p. 1850).
The full Markov chain must be irreducible.

## Usage

``` r
vcov_rel_sl(x, time, up)
```

## Arguments

- x:

  An object of class `MarkovChain`.

- time:

  A single nonnegative integer giving the discrete time point. Values
  must be less than `.Machine$integer.max`; missing and infinite values
  are not allowed.

- up:

  Integer-valued numeric vector of state positions or character vector
  of state names giving the operational states. Missing values are not
  allowed. Duplicates are ignored; empty and full state sets are valid.
  States are used in their order in `x`.

## Value

A numeric scalar containing the variance in the \\\sqrt{n}\\ limit of
the reliability estimator. Divide by the number of observed transitions
\\n\\ to approximate the variance of the estimator itself. This is a
pointwise variance, not a joint covariance across times.

## Details

The covariance of the entries of \\Q^t\\, with \\Q=P\_{UU}\\, is
computed using the component derivatives \$\$\frac{\partial
(Q^t)\_{ij}}{\partial Q\_{uv}} = \sum\_{h=0}^{t-1}
(Q^h)\_{iu}(Q^{t-1-h})\_{vj}.\$\$ Theorem 6 then sums the diagonal
variances weighted by \\a_i^2\\ and twice the covariances weighted by
\\a_i a_l\\. Every unordered pair of entries is included, interpreting
the theorem as the variance of the weighted sum; the printed
simultaneous coordinate-wise pair bounds would omit crossed pairs.
Powers and derivatives refer to the restricted block, not the
corresponding entries of the full \\P^t\\. This reference uses a dense
component covariance and is intended for small operational blocks. The
initial distribution is fixed; time zero is extended to have variance
zero.

## References

Sadek, A., and Limnios, N. (2002). *Asymptotic properties for maximum
likelihood estimators for reliability and failure rates of Markov
chains*. Communications in Statistics - Theory and Methods, 31(10),
1837–1861.
[doi:10.1081/STA-120014916](https://doi.org/10.1081/STA-120014916) .

## See also

[`rel()`](rel.md), [`vcov_rel()`](vcov_rel.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
x <- MarkovChain(c(0.4, 0.6), P)
rel(x, time = 0:3, up = 1L)
#> [1] 0.4000 0.3200 0.2560 0.2048
vcov_rel(x, time = 3L, up = 1L)
#> [1] 0.1572864
vcov_rel_sl(x, time = 3L, up = 1L)
#> [1] 0.1572864
```
