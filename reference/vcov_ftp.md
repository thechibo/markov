# Asymptotic Variance of the Failure-Time Probability

Compute the pointwise asymptotic variance of the plug-in failure-time
probability estimator for an irreducible finite-state Markov chain,
treating the initial distribution as fixed.

## Usage

``` r
vcov_ftp(x, time, up)
```

## Arguments

- x:

  An object of class `MarkovChain`.

- time:

  A single positive integer less than `.Machine$integer.max`. Missing
  and infinite values are not allowed.

- up:

  Integer-valued numeric vector of state positions or character vector
  of state names giving the operational states. Missing values are not
  allowed. Duplicates are ignored; empty and full state sets are valid.
  States are used in their order in `x`.

## Value

A numeric scalar containing the variance in the \\\sqrt{n}\\ limit of
the failure-time probability estimator. Divide by the number of observed
transitions \\n\\ to approximate the variance of the estimator itself.
This is a pointwise variance, not a joint covariance across times.

## Details

Equation (53) in Proposition 10 of Gavrilopoulos, Trevezas, and Votsi
(2026) gives the limiting variable by applying the differential
\$\$L_t(H) = a_U\left\\ \sum\_{i=0}^{t-2} Q^i H\_{UU} Q^{t-2-i} b +
Q^{t-1} H\_{UD}\mathbf{1}\_{\|D\|}\right\\,\$\$ where \\Q=P\_{UU}\\,
\\b=P\_{UD}\mathbf{1}\_{\|D\|}\\, and the sum is empty at \\t=1\\. Both
terms use the same transition estimation error; their cross-covariance
is included.

If \\G\\ contains the entry-wise coefficients of \\L_t\\, the variance
is \\g\Sigma\_{\mathbf p}g^{\mathsf T}\\, with
\\g=\operatorname{vec}(G^{\mathsf T})^{\mathsf T}\\. The implementation
contracts the full-chain multinomial row covariances using the
equivalent expression \$\$\sum\_{u\in U}\frac{1}{\pi_u}
\sum\_{j\<l}p\_{uj}p\_{ul}(G\_{uj}-G\_{ul})^2.\$\$ This includes
cross-covariances between the operational and exit blocks, avoids
constructing a dense transition covariance, and avoids subtracting
nearly equal quadratic terms. The stationary probabilities \\\pi_u\\
refer to the full chain. No new matrix dependencies are required.

## References

Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026). *Matrix
asymptotic calculus for plug-in maximum likelihood estimators in finite
Markov chains*. arXiv preprint arXiv:2607.17161.
<https://arxiv.org/abs/2607.17161>.

## See also

[`ftp()`](ftp.md), [`vcov_rel()`](vcov_rel.md), [`vcov_p()`](vcov_p.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
x <- MarkovChain(c(0.4, 0.6), P)
ftp(x, time = 1:3, up = 1L)
#> [1] 0.0800 0.0640 0.0512
vcov_ftp(x, time = 3L, up = 1L)
#> [1] 0.004369067
```
