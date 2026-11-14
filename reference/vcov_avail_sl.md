# Component-wise Reference for Availability Variance

Computes the pointwise asymptotic variance using the component-wise
covariance expansion of Sadek and Limnios (2002).

## Usage

``` r
vcov_avail_sl(x, time, up)
```

## Arguments

- x:

  An object of class `MarkovChain`.

- time:

  A single positive integer giving the time at which the asymptotic
  variance is evaluated.

- up:

  Integer-valued numeric vector of state positions or character vector
  of state names giving the operational states. Missing values are not
  allowed. Duplicates are ignored; empty and full state sets are valid.

## Value

A numeric scalar containing the pointwise asymptotic variance.

## Details

Theorem 5 (p. 1850) takes the variance of the sum of \\a_i f\_{ij}^{t}\\
over all initial states \\i\\ and operational destination states \\j\\.
The implementation adds the weighted diagonal variances and twice every
unordered pairwise covariance. The covariance entries are computed by
[`vcov_pk_sl()`](vcov_pk_sl.md) using Theorem 4 (p. 1847). Unlike
[`vcov_avail()`](vcov_avail.md), this reference does not construct a
Kronecker derivative. Arbitrary operational-state labels replace the
paper's convention that the first \\r\\ states are operational.

## References

Sadek, A., and Limnios, N. (2002). *Asymptotic properties for maximum
likelihood estimators for reliability and failure rates of Markov
chains*. Communications in Statistics - Theory and Methods, 31(10),
1837–1861.
[doi:10.1081/STA-120014916](https://doi.org/10.1081/STA-120014916) . See
Theorem 5 (p. 1850) for availability and Theorem 4 (p. 1847) for the
covariance of transition powers.
