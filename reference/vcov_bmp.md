# Asymptotic Variance of the BMP Failure Rate

Compute the pointwise asymptotic variance of the plug-in BMP
failure-rate estimator for an irreducible finite-state Markov chain with
fixed initial distribution.

## Usage

``` r
vcov_bmp(x, time, up)
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
the BMP estimator. Divide by the number of observed transitions \\n\\ to
approximate the estimator variance. This is a pointwise variance, not a
joint covariance across times.

## Details

Theorem 7 of Sadek and Limnios (2002, p. 1851) gives
\\\sigma\_\lambda^2(t)=g_t\Sigma_U g_t^{\mathsf T}\\, where \\\Sigma_U\\
is the full-chain transition covariance restricted to
operational-to-operational coordinates in row-wise order. Equation
(3.11) gives the derivative \$\$g_t =
-\frac{R(t-1)\\dR(t)-R(t)\\dR(t-1)}{R(t-1)^2}.\$\$ Both reliability
derivatives use the same transition matrix; their cross-covariance is
included. The implementation differentiates the normalized survival
recursion, an algebraically equivalent quotient rule that avoids
underflow of the reliability probabilities.

With \\G\\ the row-wise reshape of \\g_t\\, \\Q=P\_{UU}\\,
\\b_i=\sum\_{j\notin U}p\_{ij}\\, and full-chain stationary distribution
\\\pi\\, the covariance is contracted as \$\$\sum\_{i\in
U}\frac{1}{\pi_i}\left\\ \sum\_{j\<l;\\ j,l\in
U}Q\_{ij}Q\_{il}(G\_{ij}-G\_{il})^2 + b_i\sum\_{j\in
U}Q\_{ij}G\_{ij}^2\right\\.\$\$ This avoids a dense transition
covariance and subtractive cancellation.

At positive times Theorem 7 requires \\R(t-1)\>0\\; an error is raised
when this condition fails, even though [`bmpf()`](bmpf.md) defines the
rate as zero there. At time zero the variance is zero because the
initial distribution is fixed. Full operational sets also have zero
variance.

## References

Sadek, A., and Limnios, N. (2002). *Asymptotic properties for maximum
likelihood estimators for reliability and failure rates of Markov
chains*. Communications in Statistics - Theory and Methods, 31(10),
1837–1861.
[doi:10.1081/STA-120014916](https://doi.org/10.1081/STA-120014916) .

## See also

[`bmpf()`](bmpf.md), [`vcov_rel()`](vcov_rel.md),
[`vcov_p()`](vcov_p.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
x <- MarkovChain(c(0.4, 0.6), P)
bmpf(x, time = 0:3, up = 1L)
#> [1] 0.6 0.2 0.2 0.2
vcov_bmp(x, time = 3L, up = 1L)
#> [1] 0.2666667
```
