# Down-State Survival Curve of a Markov Chain

Compute the probability that a finite-state Markov chain remains in the
non-operational states from time zero through the requested time.

## Usage

``` r
dssc(x, time, down)
```

## Arguments

- x:

  An object of class `MarkovChain`.

- time:

  A nonempty nonnegative integer vector of discrete time points. Values
  must be less than `.Machine$integer.max`; missing and infinite values
  are not allowed.

- down:

  Integer-valued numeric vector of state positions or character vector
  of state names giving the non-operational states. Missing values are
  not allowed. Duplicates are ignored; empty and full state sets are
  valid. States are used in their order in `x`.

## Value

A numeric vector of the same length and order as `time` containing the
down-state survival probabilities.

## Details

For the set of down states \\D\\, the down-state survival curve is
\$\$Q(t) = \Pr\\X_0,\ldots,X_t \in D\\ = a_D P\_{DD}^t \mathbf{1}\_d,
\qquad d = \|D\|.\$\$ This is the reliability formula in Equation (2.5)
of Sadek and Limnios (2002), with the selected states changed from
operational to non-operational. The function calls [`rel()`](rel.md)
with `up = down`. Neither the restricted initial distribution nor the
transition block is renormalized. Thus \\Q(0) = \sum\_{i \in D} a_i\\;
to obtain survival conditional on starting down, supply an initial
distribution supported on `down`. Paths that leave `down` and later
return do not contribute. The curve is not generally the complement of
reliability on the up states.

## References

Sadek, A., and Limnios, N. (2002). *Asymptotic properties for maximum
likelihood estimators for reliability and failure rates of Markov
chains*. Communications in Statistics - Theory and Methods, 31(10),
1837–1861.
[doi:10.1081/STA-120014916](https://doi.org/10.1081/STA-120014916) .

## See also

[`rel()`](rel.md), [`mttr()`](mttr.md), [`vcov_dssc()`](vcov_dssc.md),
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
