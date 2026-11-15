# Asymptotic Confidence Intervals for Mean Time to Failure

Compute two-sided, pointwise normal-approximation confidence intervals
using [`mttf()`](mttf.md) and [`var_mttf()`](var_mttf.md).

## Usage

``` r
ci_mttf(x, up, n, level = 0.95)
```

## Arguments

- x:

  An object of class `MarkovChain`.

- up:

  Nonempty integer-valued numeric vector of state positions, or
  character vector of state names, giving the operational states.
  Entries must be distinct and select a proper subset of the states.
  Input order is ignored; states are used in their order in `x`.

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
confidence limits. There is one row.

## Details

The asymptotic normality result is Proposition 11 of Gavrilopoulos,
Trevezas, and Votsi (2026). The domain restrictions of
[`var_mttf()`](var_mttf.md) apply.

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

[`mttf()`](mttf.md), [`var_mttf()`](var_mttf.md),
[`ci_rel()`](ci_rel.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE)
x <- MarkovChain(c(0.4, 0.6), P)
ci_mttf(x, up = 1L, n = 1000)
#>   estimate std.error    lower    upper
#> 1        2 0.1632993 1.679939 2.320061
```
