# Failure-Time Probability of a Markov Chain

Compute the probability that a finite-state Markov chain first enters a
non-operational state at the requested positive time.

## Usage

``` r
ftp(x, time, up)
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
failure-time probabilities.

## Details

Let \\U\\ be the operational states, \\D\\ their complement, and \\T_D =
\inf\\t \geq 0 : X_t \in D\\\\. Equation (44) of Gavrilopoulos,
Trevezas, and Votsi (2026) gives \$\$f(t) = \Pr(T_D=t) = a_U Q^{t-1} b,
\qquad Q=P\_{UU}, \quad b=P\_{UD}\mathbf{1}\_{\|D\|}, \quad t\geq1.\$\$
Neither the initial distribution nor the restricted block is
renormalized. Paths starting down have failure time zero and do not
contribute at the positive times supported here. Their probability is
\\\sum\_{i\in D}a_i\\. The exit vector is summed directly from
\\P\_{UD}\\, avoiding subtraction of nearly equal reliability
probabilities. Powers are computed by [`Ppowers()`](Ppowers.md). Empty
and full operational sets give zero at positive times.

## References

Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026). *Matrix
asymptotic calculus for plug-in maximum likelihood estimators in finite
Markov chains*. arXiv preprint arXiv:2607.17161.
<https://arxiv.org/abs/2607.17161>.

## See also

[`rel()`](rel.md), [`vcov_ftp()`](vcov_ftp.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
x <- MarkovChain(c(0.4, 0.6), P)
ftp(x, time = 1:3, up = 1L)
#> [1] 0.0800 0.0640 0.0512
vcov_ftp(x, time = 3L, up = 1L)
#> [1] 0.004369067
```
