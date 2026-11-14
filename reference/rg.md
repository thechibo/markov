# RG Failure Rate of a Markov Chain

Compute the discrete-time RG failure rate as the logarithmic decrease in
reliability over one time step.

## Usage

``` r
rg(x, time, up)
```

## Arguments

- x:

  An object of class `MarkovChain`.

- time:

  A nonempty nonnegative integer vector of discrete time points. Values
  must be less than `.Machine$integer.max`; missing and infinite values
  are not allowed.

- up:

  Integer-valued numeric vector of state positions or character vector
  of state names giving the operational states. Missing values are not
  allowed. Duplicates are ignored; empty and full state sets are valid.
  States are used in their order in `x`.

## Value

A numeric vector of the same length and order as `time`, containing RG
rates, with `Inf` and `NaN` as described above.

## Details

Equation (2.9) of Sadek and Limnios (2002, p. 1842) defines
\$\$r(t)=\log\\R(t-1)/R(t)\\, \quad t\geq1, \qquad r(0)=-\log R(0).\$\$
When the reliability ratio is defined, this equals
\\-\log(1-\lambda(t))\\, with BMP rate \\\lambda(t)\\. The
implementation reuses the normalized survival recursion underlying
[`bmp()`](bmp.md), avoiding underflow of the raw reliability powers. It
uses [`log1p()`](https://rdrr.io/r/base/Log.html) for small failure
probabilities and the directly computed conditional survival probability
when the BMP rate is close to one. The initial distribution is treated
as fixed.

If \\R(t-1)\>0\\ but \\R(t)=0\\, the rate is `Inf`. If both are zero,
the defining ratio is undefined and `NaN` is returned, rather than
applying the zero-rate convention of [`bmp()`](bmp.md). At time zero,
\\R(0)=0\\ gives `Inf`. Full operational sets have rate zero.

## References

Sadek, A., and Limnios, N. (2002). *Asymptotic properties for maximum
likelihood estimators for reliability and failure rates of Markov
chains*. Communications in Statistics - Theory and Methods, 31(10),
1837–1861.
[doi:10.1081/STA-120014916](https://doi.org/10.1081/STA-120014916) .

## See also

[`bmp()`](bmp.md), [`rel()`](rel.md), [`vcov_rg()`](vcov_rg.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
x <- MarkovChain(c(0.4, 0.6), P)
rg(x, time = 0:3, up = 1L)
#> [1] 0.9162907 0.2231436 0.2231436 0.2231436
vcov_rg(x, time = 3L, up = 1L)
#> [1] 0.4166667
```
