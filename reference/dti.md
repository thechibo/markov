# Discrete-Time Intensity of Hitting Down

Compute the probability of an operational-to-down transition at a given
positive time in a finite-state Markov chain.

## Usage

``` r
dti(x, time, up)
```

## Arguments

- x:

  An object of class `MarkovChain`.

- time:

  A nonempty positive integer vector of discrete time points, each less
  than `.Machine$integer.max`. Missing and infinite values are not
  allowed.

- up:

  Integer-valued numeric vector of state positions or character vector
  of state names giving the operational states. Missing values are not
  allowed. Duplicates are ignored; empty and full state sets are valid.
  States are used in their order in `x`.

## Value

A numeric vector of the same length and order as `time` containing the
discrete-time intensities.

## Details

Equation (48) of Gavrilopoulos, Trevezas, and Votsi (2026) gives
\$\$r(t) = \Pr(X\_{t-1}\in U, X_t\in D) = a P^{t-1} J_U P e_D, \qquad
t\geq1,\$\$ where \\U\\ is the operational set, \\D\\ its complement,
\\J_U\\ is the diagonal operational-state indicator matrix, and \\e_D\\
is the down-state indicator vector. The full transition matrix evolves
the initial distribution to time \\t-1\\; paths may have failed and
recovered before that time. This differs from the first-entry
probability returned by [`ftp()`](ftp.md). Paths starting down can
contribute after recovery. The initial distribution is fixed and is not
renormalized.

The vector \\J_U P e_D\\ is computed by summing down-transition
probabilities directly and setting its down-state entries to zero.
Powers are computed once by [`Ppowers()`](Ppowers.md). Empty and full
operational sets give zero intensity.

## References

Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026). *Matrix
asymptotic calculus for plug-in maximum likelihood estimators in finite
Markov chains*. arXiv preprint arXiv:2607.17161.
<https://arxiv.org/abs/2607.17161>.

## See also

[`ftp()`](ftp.md), [`avail()`](avail.md), [`vcov_dti()`](vcov_dti.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
x <- MarkovChain(c(0.4, 0.6), P)
dti(x, time = 1:3, up = 1L)
#> [1] 0.08 0.10 0.11
vcov_dti(x, time = 3L, up = 1L)
#> [1] 0.05826667
```
