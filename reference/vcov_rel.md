# Asymptotic Variance of Reliability

Compute the pointwise asymptotic variance of the plug-in reliability
estimator for an irreducible finite-state Markov chain, treating the
initial distribution as fixed.

## Usage

``` r
vcov_rel(x, time, up)
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
the reliability estimator. Divide by the number of observed transitions
\\n\\ to approximate the variance of the estimator itself. This is a
pointwise variance, not a joint covariance across times.

## Details

Proposition 10, Equation (52), of Gavrilopoulos, Trevezas, and Votsi
(2026) gives the variance as \\D_t \Sigma_U D_t^{\mathsf T}\\, where
\$\$D_t = \sum\_{i=0}^{t-1} (a_U Q^i) \otimes
(Q^{t-1-i}\mathbf{1}\_r)^{\mathsf T}.\$\$ Here \\Q=P\_{UU}\\ and
\\\Sigma_U\\ selects the operational-transition coordinates from the
full-chain transition covariance in row-wise order. The Kronecker
operator is contracted with the initial distribution and vector of ones
before forming the variance, avoiding a dense \\r^2 \times r^2\\ power
derivative. The full-chain stationary distribution is used; the
operational block is not renormalized. At time zero the derivative and
variance are zero.

## References

Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026). *Matrix
asymptotic calculus for plug-in maximum likelihood estimators in finite
Markov chains*. arXiv preprint arXiv:2607.17161.
<https://arxiv.org/abs/2607.17161>.

## See also

[`rel()`](rel.md), [`vcov_rel_sl()`](vcov_rel_sl.md),
[`vcov_p()`](vcov_p.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
x <- MarkovChain(c(0.4, 0.6), P)
rel(x, time = 0:3, up = 1L)
#> [1] 0.4000 0.3200 0.2560 0.2048
vcov_rel(x, time = 3L, up = 1L)
#> [1] 0.1572864
vcov_rel_sl(x, time = 3L, up = 1L)
#> [1] 0.1572864
```
