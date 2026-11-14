# Maintainability of a Markov Chain

Compute the probability that a finite-state Markov chain has reached an
operational state by the requested time, including time zero.

## Usage

``` r
main(x, time, down)
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
maintainability probabilities.

## Details

Let \\D\\ be the set of non-operational states and \\T_U = \inf\\t \geq
0 : X_t \notin D\\\\. Maintainability is \$\$M(t) = \Pr(T_U \leq t) =
1 - Q(t),\$\$ where \\Q(t)\\ is the down-state survival curve. This
function computes `1 - dssc(x, time, down)`, as in Proposition 10 of
Gavrilopoulos, Trevezas, and Votsi (2026).

The initial distribution is not renormalized. Consequently, \\M(0) = 1 -
\sum\_{i \in D} a_i\\: paths starting operational have recovery time
zero. To obtain maintainability conditional on starting down, supply an
initial distribution supported on `down`. Paths that recover and
subsequently fail still contribute to maintainability.

## References

Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026). *Matrix
asymptotic calculus for plug-in maximum likelihood estimators in finite
Markov chains*. arXiv preprint arXiv:2607.17161.
<https://arxiv.org/abs/2607.17161>.

## See also

[`dssc()`](dssc.md), [`mttr()`](mttr.md), [`vcov_main()`](vcov_main.md),
[`vcov_main_sl()`](vcov_main_sl.md)

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
