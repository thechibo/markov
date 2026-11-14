# Asymptotic Variance of Maintainability

Compute the pointwise asymptotic variance of the plug-in maintainability
estimator for an irreducible finite-state Markov chain, treating the
initial distribution as fixed.

## Usage

``` r
vcov_main(x, time, down)
```

## Arguments

- x:

  An object of class `MarkovChain`.

- time:

  A single nonnegative integer giving the discrete time point. Values
  must be less than `.Machine$integer.max`; missing and infinite values
  are not allowed.

- down:

  Integer-valued numeric vector of state positions or character vector
  of state names giving the non-operational states. Missing values are
  not allowed. Duplicates are ignored; empty and full state sets are
  valid. States are used in their order in `x`.

## Value

A numeric scalar containing the variance in the \\\sqrt{n}\\ limit of
the maintainability estimator. Divide by the number of observed
transitions \\n\\ to approximate the variance of the estimator itself.
This is a pointwise variance, not a joint covariance across times.

## Details

Since \\\widehat M(t) - M(t) = -(\widehat Q(t) - Q(t))\\, the asymptotic
variance of maintainability equals that of down-state survival.
Proposition 10 of Gavrilopoulos, Trevezas, and Votsi (2026) states this
sign change after Equation (54). The variance is obtained from Equation
(52) with the down-state block substituted for the operational block.
This function calls [`vcov_dssc()`](vcov_dssc.md) without changing its
result. The full-chain stationary distribution is used and the variance
at time zero is zero.

## References

Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026). *Matrix
asymptotic calculus for plug-in maximum likelihood estimators in finite
Markov chains*. arXiv preprint arXiv:2607.17161.
<https://arxiv.org/abs/2607.17161>.

## See also

[`main()`](main.md), [`vcov_dssc()`](vcov_dssc.md),
[`vcov_main_sl()`](vcov_main_sl.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
x <- MarkovChain(c(0, 1), P)
main(x, time = 0:3, down = 2L)
#> [1] 0.000 0.300 0.510 0.657
vcov_main(x, time = 3L, down = 2L)
#> [1] 1.134472
vcov_main_sl(x, time = 3L, down = 2L)
#> [1] 1.134472
```
