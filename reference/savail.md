# Stationary Availability of a Markov Chain

Compute the stationary probability of being operational for an
irreducible finite-state Markov chain.

## Usage

``` r
savail(x, up)
```

## Arguments

- x:

  An object of class `MarkovChain`.

- up:

  Integer-valued numeric vector of state positions or character vector
  of state names giving the operational states. Missing values are not
  allowed. Duplicates are ignored; empty and full state sets are valid.
  States are used in their order in `x`.

## Value

A numeric scalar containing the stationary availability.

## Details

Equation (49) of Gavrilopoulos, Trevezas, and Votsi (2026) gives
\$\$A\_\infty = \pi e_U = \sum\_{i\in U}\pi_i,\$\$ where \\\pi\\ is the
stationary distribution of the full transition matrix and \\e_U\\ is the
operational-state indicator vector. This function uses the same
stationary solver as [`stat()`](stat.md). Both this function and
[`vcov_savail()`](vcov_savail.md) require an irreducible chain, so the
stationary distribution is unique. The result does not depend on the
initial distribution, although the complete `MarkovChain` is validated.

Periodic chains are supported. For them, stationary availability need
not equal the limit of [`avail()`](avail.md) from an arbitrary initial
distribution. Empty and full operational sets give zero and one,
respectively.

## References

Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026). *Matrix
asymptotic calculus for plug-in maximum likelihood estimators in finite
Markov chains*. arXiv preprint arXiv:2607.17161.
<https://arxiv.org/abs/2607.17161>.

## See also

[`stat()`](stat.md), [`avail()`](avail.md),
[`vcov_savail()`](vcov_savail.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
x <- MarkovChain(c(0.4, 0.6), P)
savail(x, up = 1L)
#> [1] 0.6
vcov_savail(x, up = 1L)
#> [1] 0.72
```
