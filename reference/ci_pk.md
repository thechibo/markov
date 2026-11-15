# Confidence Intervals for k-Step Transition Probabilities

Compute a two-sided, pointwise normal-approximation confidence interval
for each entry of the k-step transition matrix \\P^k\\.

## Usage

``` r
ci_pk(x, k, n, level = 0.95)
```

## Arguments

- x:

  A [`MarkovChain()`](MarkovChain.md) or [`StochMat()`](StochMat.md)
  containing the transition matrix estimated from a single observed
  trajectory. The estimated chain must be irreducible.

- k:

  A single positive integer giving the number of steps.

- n:

  A single finite positive integer giving the number of observed
  one-step transitions used to estimate the transition matrix.

- level:

  A single finite numeric confidence level strictly between zero and
  one. Defaults to `0.95`; missing values are not allowed.

## Value

A data frame with one row per transition, ordered by origin state and
then by destination state in matrix order. Character columns `from` and
`to` identify the states. Numeric columns `estimate`, `std.error`,
`lower`, and `upper` give the plug-in probability, its estimated
standard error, and its confidence limits.

## Details

The estimates are the entries of the fitted transition matrix raised to
the kth power. The intervals use the diagonal of
[`vcov_pk()`](vcov_pk.md), based on Proposition 4, Equation (26), of
Gavrilopoulos, Trevezas, and Votsi (2026). The normal approximation
follows Remark 3, as in [`ci_p()`](ci_p.md). For `k = 1`, the two
functions give the same intervals.

## Pointwise coverage

The confidence level applies to each entry separately. Bounds may extend
outside \\\[0,1\]\\. A zero estimated variance gives a collapsed
interval. Nominal asymptotic coverage requires a positive limiting
variance.

## References

Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026). *Matrix
asymptotic calculus for plug-in maximum likelihood estimators in finite
Markov chains*. arXiv preprint arXiv:2607.17161v2.
<https://arxiv.org/abs/2607.17161v2>.

## See also

[`ci_p()`](ci_p.md), [`vcov_pk()`](vcov_pk.md),
[`Ppowers()`](Ppowers.md)

## Examples

``` r
P <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE))
ci_pk(P, k = 3, n = 1000)
#>   from to estimate  std.error     lower     upper
#> 1    1  1    0.650 0.02387467 0.6032065 0.6967935
#> 2    1  2    0.350 0.02387467 0.3032065 0.3967935
#> 3    2  1    0.525 0.02811250 0.4699005 0.5800995
#> 4    2  2    0.475 0.02811250 0.4199005 0.5300995
```
