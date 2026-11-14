# Availability of a Markov Chain

Compute the discrete-time availability function of a system whose
degradation/repair dynamics are modeled by a finite-state Markov chain.
Availability at time \\t\\ is defined as \$\$ A(t) = \Pr\\ X_t \in
\mathcal{U} \\ \$\$ where \\X_t\\ is the state of the chain at time
\\t\\ and \\\mathcal{U}\\ is the set of *up* (operational) states.

## Usage

``` r
avail(x, time, up)
```

## Arguments

- x:

  An object of class `MarkovChain`.

- time:

  A nonempty positive integer vector giving the discrete time points
  \\t\\ at which availability \\A(t)\\ is evaluated.

- up:

  Integer-valued numeric vector of state positions or character vector
  of state names giving the operational states. Missing values are not
  allowed. Duplicates are ignored; empty and full state sets are valid.

## Value

A numeric vector of the same length as `time`, where the \\j\\-th
element is the probability that the system is available (in an up state)
at `time[j]`.

## Details

For each time \\t\\, the function computes \$\$ A(t) = a P^t
\mathbf{1}\_{\mathcal{U}}, \$\$ where \\a\\ is the initial distribution,
\\P\\ is the transition matrix, and \\\mathbf{1}\_{\mathcal{U}}\\ is the
indicator vector of the operational states.

Duplicate values in `up` are ignored. Powers are computed by repeated
real matrix multiplication with [`Ppowers()`](Ppowers.md).

## See also

[`MarkovChain()`](MarkovChain.md), [`Ppowers()`](Ppowers.md),
[`vcov_avail()`](vcov_avail.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
x <- MarkovChain(c(0.4, 0.6), P)
avail(x, time = 1:3, up = 1L)
#> [1] 0.500 0.550 0.575
vcov_avail(x, time = 3L, up = 1L)
#> [1] 0.5953125
```
