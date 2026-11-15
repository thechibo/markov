# Confidence Intervals for Transition Probabilities

Compute a two-sided, pointwise normal-approximation confidence interval
for each entry of the one-step transition matrix.

## Usage

``` r
ci_p(x, n, level = 0.95)
```

## Arguments

- x:

  A [`MarkovChain()`](MarkovChain.md) or [`StochMat()`](StochMat.md)
  containing the transition matrix estimated from a single observed
  trajectory. The estimated chain must be irreducible.

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

The intervals use the diagonal of the row-wise asymptotic covariance
matrix from [`vcov_p()`](vcov_p.md), as in Theorem 1, Equation (2), of
Gavrilopoulos, Trevezas, and Votsi (2026). The standard errors divide
the square roots of these variances by \\\sqrt{n}\\ before applying the
normal quantile in Remark 3.

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

[`vcov_p()`](vcov_p.md), [`ci_pk()`](ci_pk.md)

## Examples

``` r
P <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE))
ci_p(P, n = 1000)
#>   from to estimate  std.error     lower     upper
#> 1    1  1      0.8 0.01632993 0.7679939 0.8320061
#> 2    1  2      0.2 0.01632993 0.1679939 0.2320061
#> 3    2  1      0.3 0.02291288 0.2550916 0.3449084
#> 4    2  2      0.7 0.02291288 0.6550916 0.7449084
```
