# Component-wise Reference for Down-State Survival Variance

Compute the pointwise asymptotic variance of down-state survival using
the component-wise covariance expansion of Sadek and Limnios (2002). The
full Markov chain must be irreducible.

## Usage

``` r
vcov_dssc_sl(x, time, down)
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

This function calls [`vcov_rel_sl()`](vcov_rel_sl.md) with `up = down`,
applying Theorem 6 (p. 1850) to the down-state block. It sums the
weighted diagonal variances and twice every unordered pairwise
covariance of the entries of \\P\_{DD}^t\\. As in
[`vcov_rel_sl()`](vcov_rel_sl.md), all unordered pairs are included to
implement the variance of the weighted sum; the printed simultaneous
coordinate-wise pair bounds would omit crossed pairs. The initial
distribution is fixed and its restriction to `down` is not renormalized.
Time zero has variance zero. This reference constructs a dense component
covariance and is intended for small down-state blocks.

## References

Sadek, A., and Limnios, N. (2002). *Asymptotic properties for maximum
likelihood estimators for reliability and failure rates of Markov
chains*. Communications in Statistics - Theory and Methods, 31(10),
1837–1861.
[doi:10.1081/STA-120014916](https://doi.org/10.1081/STA-120014916) .

## See also

[`dssc()`](dssc.md), [`vcov_dssc()`](vcov_dssc.md),
[`vcov_rel_sl()`](vcov_rel_sl.md)

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
