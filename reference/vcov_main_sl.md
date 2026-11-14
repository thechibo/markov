# Component-wise Reference for Maintainability Variance

Compute the pointwise asymptotic variance of maintainability using the
component-wise covariance expansion of Sadek and Limnios (2002). The
full Markov chain must be irreducible.

## Usage

``` r
vcov_main_sl(x, time, down)
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
the maintainability estimator. Divide by the number of observed
transitions \\n\\ to approximate the variance of the estimator itself.
This is a pointwise variance, not a joint covariance across times.

## Details

This function calls [`vcov_dssc_sl()`](vcov_dssc_sl.md) because
subtracting down-state survival from one leaves its variance unchanged.
The underlying reference applies Theorem 6 (p. 1850) of Sadek and
Limnios (2002) to the down-state block, including every unordered
pairwise covariance as documented in [`vcov_rel_sl()`](vcov_rel_sl.md).
The initial distribution is fixed and is not renormalized on `down`.
Time zero has variance zero. The reference constructs a dense component
covariance and is intended for small down-state blocks.

## References

Sadek, A., and Limnios, N. (2002). *Asymptotic properties for maximum
likelihood estimators for reliability and failure rates of Markov
chains*. Communications in Statistics - Theory and Methods, 31(10),
1837–1861.
[doi:10.1081/STA-120014916](https://doi.org/10.1081/STA-120014916) .

## See also

[`main()`](main.md), [`vcov_main()`](vcov_main.md),
[`vcov_dssc_sl()`](vcov_dssc_sl.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
x <- MarkovChain(c(0, 1), P)
main(x, time = 0:3, down = 2L)
#> [1] 0.000 0.300 0.510 0.657
vcov_main(x, time = 3L, down = 2L)
#> [1] 1.134472
vcov_main_sl(x, time = 3L, down = 2L)
#> [1] 1.134472
```
