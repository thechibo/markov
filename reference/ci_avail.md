# Asymptotic Confidence Intervals for Availability

Compute two-sided, pointwise normal-approximation confidence intervals
using [`avail()`](avail.md) and [`vcov_avail()`](vcov_avail.md).

## Usage

``` r
ci_avail(x, time, up, n, level = 0.95)
```

## Arguments

- x:

  An object of class `MarkovChain`.

- time:

  A nonempty positive integer vector giving the discrete time points
  \\t\\ at which availability \\A(t)\\ is evaluated.

- up:

  Integer-valued numeric vector of state positions or character vector
  of state names giving the operational states. Missing values are not
  allowed. Duplicates are ignored; empty and full state sets are valid.

- n:

  A single finite positive integer: the number of observed transitions
  used to estimate `x$P` (a trajectory of length `n + 1`). Missing
  values are not allowed.

- level:

  A single finite numeric confidence level strictly between zero and
  one. Defaults to `0.95`; missing values are not allowed.

## Value

A data frame with numeric columns `estimate`, `std.error`, `lower`, and
`upper`, giving the point estimate, estimated standard error, and
confidence limits. There is one row per requested time, preserving order
and duplicates.

## Details

The asymptotic normality result is Proposition 9 of Gavrilopoulos,
Trevezas, and Votsi (2026). The domain restrictions of
[`vcov_avail()`](vcov_avail.md) apply.

## Normal approximation

The interval uses the first-order normal approximation for a scalar
plug-in estimator described in Remark 3 of Gavrilopoulos, Trevezas, and
Votsi (2026). Supply a chain whose transition matrix is the plug-in
maximum likelihood estimate from the observed trajectory. The chain must
be irreducible. The initial distribution is treated as fixed, so its
estimation uncertainty is not included. The variance function returns
the variance in the \\\sqrt{n}\\ limit, giving the interval
\$\$\widehat{\theta} \mathbin{\pm}
z\_{(1+\mathrm{level})/2}\sqrt{\widehat{\sigma}^2/n}.\$\$ Bounds are not
clipped to the parameter space. Intervals at multiple times are
pointwise, not a simultaneous confidence band. Nominal asymptotic
coverage requires a positive limiting variance and the assumptions of
the linked variance function. Zero estimated variance produces a
collapsed interval; this does not establish nominal coverage at a
degenerate boundary. Nonfinite estimates or variances and negative
variances produce an error.

## References

Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026). *Matrix
asymptotic calculus for plug-in maximum likelihood estimators in finite
Markov chains*. arXiv preprint arXiv:2607.17161v2.
<https://arxiv.org/abs/2607.17161v2>.

## See also

[`avail()`](avail.md), [`vcov_avail()`](vcov_avail.md),
[`ci_rel()`](ci_rel.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE)
x <- MarkovChain(c(0.4, 0.6), P)
ci_avail(x, time = 1:3, up = 1L, n = 1000)
#>   estimate  std.error     lower     upper
#> 1    0.500 0.01522060 0.4701682 0.5298318
#> 2    0.550 0.02160247 0.5076599 0.5923401
#> 3    0.575 0.02439903 0.5271788 0.6228212
```
