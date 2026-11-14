# Reliability of a Markov Chain

Compute the probability that a finite-state Markov chain remains in the
operational states from time zero through the requested time.

## Usage

``` r
rel(x, time, up)
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

A numeric vector of the same length and order as `time` containing the
reliability probabilities.

## Details

Equation (2.5) of Sadek and Limnios (2002) gives \$\$R(t) = a_U Q^t
\mathbf{1}\_r, \qquad Q = P\_{UU},\$\$ where \\a_U\\ is the initial
distribution restricted to the operational states and \\r\\ is their
number. Neither \\a_U\\ nor \\Q\\ is renormalized. In particular, \\R(0)
= \sum\_{i \in U} a_i\\. Paths that leave `up` and subsequently return
do not contribute to reliability. Powers of the restricted block are
computed by [`Ppowers()`](Ppowers.md), without an eigendecomposition.

## References

Sadek, A., and Limnios, N. (2002). *Asymptotic properties for maximum
likelihood estimators for reliability and failure rates of Markov
chains*. Communications in Statistics - Theory and Methods, 31(10),
1837–1861.
[doi:10.1081/STA-120014916](https://doi.org/10.1081/STA-120014916) .

## See also

[`avail()`](avail.md), [`mttf()`](mttf.md), [`vcov_rel()`](vcov_rel.md),
[`vcov_rel_sl()`](vcov_rel_sl.md)

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
