# Asymptotic Variance of Stationary Availability

Compute the asymptotic variance of the plug-in stationary availability
estimator for an irreducible finite-state Markov chain.

## Usage

``` r
vcov_savail(x, up)
```

## Arguments

- x:

  An object of class `MarkovChain`.

- up:

  Integer-valued numeric vector of state positions or character vector
  of state names giving the operational states. Missing values are not
  allowed. Duplicates are ignored; empty and full state sets are valid.
  States are used in their order in `x`.

## Value

A numeric scalar containing the variance in the \\\sqrt{n}\\ limit of
the stationary availability estimator. Divide by the number of observed
transitions \\n\\ to approximate the variance of the estimator itself.

## Details

Equation (60) in Proposition 12 of Gavrilopoulos, Trevezas, and Votsi
(2026) gives the limiting variable \\W\_\pi e_U\\. Consequently,
\$\$\sigma^2\_{A\_\infty} = e_U^{\mathsf T}\Sigma\_\pi e_U,\$\$ where
\\\Sigma\_\pi\\ is the stationary-distribution asymptotic covariance
from Proposition 7, computed by [`vcov_stat()`](vcov_stat.md). The
function sums all entries of its operational-state submatrix, including
the off-diagonal covariances. Since the stationary probabilities sum to
one, the complementary state set has the same variance; the smaller set
is used to reduce summation and cancellation. Empty and full sets return
exactly zero after checking irreducibility.

The result does not depend on the initial distribution. This is the
asymptotic variance of the transition-MLE plug-in estimator of
stationary availability. The numerical conditioning limitations of
[`vcov_stat()`](vcov_stat.md) also apply here.

## References

Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026). *Matrix
asymptotic calculus for plug-in maximum likelihood estimators in finite
Markov chains*. arXiv preprint arXiv:2607.17161.
<https://arxiv.org/abs/2607.17161>.

## See also

[`savail()`](savail.md), [`vcov_stat()`](vcov_stat.md),
[`vcov_avail()`](vcov_avail.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
x <- MarkovChain(c(0.4, 0.6), P)
savail(x, up = 1L)
#> [1] 0.6
vcov_savail(x, up = 1L)
#> [1] 0.72
```
