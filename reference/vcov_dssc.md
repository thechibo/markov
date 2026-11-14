# Asymptotic Variance of Down-State Survival

Compute the pointwise asymptotic variance of the plug-in down-state
survival estimator for an irreducible finite-state Markov chain,
treating the initial distribution as fixed.

## Usage

``` r
vcov_dssc(x, time, down)
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
the down-state survival estimator. Divide by the number of observed
transitions \\n\\ to approximate the variance of the estimator itself.
This is a pointwise variance, not a joint covariance across times.

## Details

Proposition 10, Equation (54), of Gavrilopoulos, Trevezas, and Votsi
(2026) gives the down-state survival limit. Its variance is Equation
(52) with the operational block replaced by the down-state block:
\$\$\sigma_Q^2(t) = D_t \Sigma_D D_t^{\mathsf T}, \qquad D_t =
\sum\_{i=0}^{t-1} (a_D P\_{DD}^i) \otimes
(P\_{DD}^{t-1-i}\mathbf{1}\_d)^{\mathsf T}.\$\$ The covariance
\\\Sigma_D\\ selects the down-to-down transition coordinates in row-wise
order from the full-chain transition covariance. The function calls
[`vcov_rel()`](vcov_rel.md) with `up = down`, retaining the full-chain
stationary distribution and the unnormalized restricted initial
distribution. The variance at time zero is zero.

## References

Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026). *Matrix
asymptotic calculus for plug-in maximum likelihood estimators in finite
Markov chains*. arXiv preprint arXiv:2607.17161.
<https://arxiv.org/abs/2607.17161>.

## See also

[`dssc()`](dssc.md), [`vcov_rel()`](vcov_rel.md),
[`vcov_dssc_sl()`](vcov_dssc_sl.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
x <- MarkovChain(c(0.4, 0.6), P)
dssc(x, time = 0:3, down = 2L)
#> [1] 0.6000 0.4200 0.2940 0.2058
vcov_dssc(x, time = 3L, down = 2L)
#> [1] 0.4084101
vcov_dssc_sl(x, time = 3L, down = 2L)
#> [1] 0.4084101
```
