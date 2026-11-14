# Variance of availability

Compute the pointwise asymptotic variance of the plug-in availability
estimator for an irreducible finite-state Markov chain, treating the
initial distribution as fixed.

## Usage

``` r
vcov_avail(x, time, up)
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

A numeric scalar containing the variance in the \\\sqrt{m}\\ limit of
the availability estimator, where \\m\\ is the sample size. Divide by
\\m\\ to approximate the variance of the estimator itself. This is a
pointwise variance, not a joint covariance across times.

## Details

This function evaluates the diagonal entry of Equation (29) for the
requested time, using the operational-state indicator specified in
Proposition 9 of Gavrilopoulos, Trevezas, and Votsi (2026). The
derivative is assembled as a sum of Kronecker products, using real
transition powers computed once by [`Ppowers()`](Ppowers.md). No
eigendecomposition or inverse eigenvector matrix is required.

## References

Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026). *Matrix
asymptotic calculus for plug-in maximum likelihood estimators in finite
Markov chains*. arXiv preprint arXiv:2607.17161.
<https://arxiv.org/abs/2607.17161>.

## See also

[`avail()`](avail.md), [`vcov_p()`](vcov_p.md),
[`vcov_pk()`](vcov_pk.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
x <- MarkovChain(c(0.4, 0.6), P)
avail(x, time = 1:3, up = 1L)
#> [1] 0.500 0.550 0.575
vcov_avail(x, time = 3L, up = 1L)
#> [1] 0.5953125
```
