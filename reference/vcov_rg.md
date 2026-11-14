# Asymptotic Variance of the RG Failure Rate

Compute the pointwise asymptotic variance of the plug-in RG failure-rate
estimator for an irreducible finite-state Markov chain with fixed
initial distribution.

## Usage

``` r
vcov_rg(x, time, up)
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
the RG estimator. Divide by the number of observed transitions \\n\\ to
approximate the estimator variance. This is a pointwise variance, not a
joint covariance across times.

## Details

Corollary 2 of Sadek and Limnios (2002, pp. 1851–1852) applies the delta
method to \\r(t)=-\log(1-\lambda(t))\\, giving
\$\$\sigma_r^2(t)=\frac{\sigma\_\lambda^2(t)}{(1-\lambda(t))^2}.\$\$
This function reuses [`vcov_bmp()`](vcov_bmp.md) for
\\\sigma\_\lambda^2(t)\\ and the normalized survival recursion for
\\1-\lambda(t)=R(t)/R(t-1)\\. The scaling is performed by two successive
divisions, avoiding premature underflow from squaring a small survival
probability.

At positive times, both preceding and current survival must be positive.
An error is raised for infinite or undefined RG rates, where Corollary 2
does not provide a finite-rate delta-method limit. At time zero the
variance is zero when \\R(0)\>0\\, because the initial distribution is
fixed. Full operational sets also have variance zero. Very small
survival ratios can produce an asymptotic variance exceeding the numeric
range.

## References

Sadek, A., and Limnios, N. (2002). *Asymptotic properties for maximum
likelihood estimators for reliability and failure rates of Markov
chains*. Communications in Statistics - Theory and Methods, 31(10),
1837–1861.
[doi:10.1081/STA-120014916](https://doi.org/10.1081/STA-120014916) .

## See also

[`rgf()`](rgf.md), [`vcov_bmp()`](vcov_bmp.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
x <- MarkovChain(c(0.4, 0.6), P)
rgf(x, time = 0:3, up = 1L)
#> [1] 0.9162907 0.2231436 0.2231436 0.2231436
vcov_rg(x, time = 3L, up = 1L)
#> [1] 0.4166667
```
