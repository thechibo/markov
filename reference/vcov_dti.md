# Asymptotic Variance of the Discrete-Time Intensity

Compute the pointwise asymptotic variance of the plug-in discrete-time
intensity estimator for an irreducible finite-state Markov chain,
treating the initial distribution as fixed.

## Usage

``` r
vcov_dti(x, time, up)
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
the intensity estimator. Divide by the number of observed transitions
\\n\\ to approximate the variance of the estimator itself. This is a
pointwise variance, not a joint covariance across times.

## Details

Equation (59) in Proposition 12 of Gavrilopoulos, Trevezas, and Votsi
(2026) applies the differential \$\$L_t(H) = a\left\\ \sum\_{i=0}^{t-2}
P^i H P^{t-2-i} J_U P + P^{t-1} J_U H\right\\e_D.\$\$ The sum is empty
at \\t=1\\. Both the distribution at time \\t-1\\ and the final
transition contribute; their cross-covariance is included. At later
times, transitions from down states can affect the variance through the
first term.

With \\G\\ the matrix of entry-wise coefficients of \\L_t\\, let
\\g=\operatorname{vec}(G^{\mathsf T})^{\mathsf T}\\ denote the row-wise
derivative. The variance \\g\Sigma\_{\mathbf p}g^{\mathsf T}\\ is
computed using the equivalent full-chain row-covariance contraction
\$\$\sum_u \frac{1}{\pi_u}\sum\_{j\<l}
p\_{uj}p\_{ul}(G\_{uj}-G\_{ul})^2,\$\$ where \\\pi\\ is the stationary
distribution of the full chain. This avoids constructing a dense
transition covariance or subtracting nearly equal quadratic terms. The
derivative is assembled using outer products, without a full
matrix-power Jacobian.

## References

Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026). *Matrix
asymptotic calculus for plug-in maximum likelihood estimators in finite
Markov chains*. arXiv preprint arXiv:2607.17161.
<https://arxiv.org/abs/2607.17161>.

## See also

[`dti()`](dti.md), [`vcov_ftp()`](vcov_ftp.md), [`vcov_p()`](vcov_p.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
x <- MarkovChain(c(0.4, 0.6), P)
dti(x, time = 1:3, up = 1L)
#> [1] 0.08 0.10 0.11
vcov_dti(x, time = 3L, up = 1L)
#> [1] 0.05826667
```
