# Simulate and Estimate Markov Chain Functionals

Simulate independent trajectories, fit a Markov chain to each
trajectory, and apply a function to each fit.

## Usage

``` r
fit_parallel(x, fun, k = 1000, m = 500, ..., n_workers = 2L)
```

## Arguments

- x:

  A `MarkovChain` object.

- fun:

  A function accepting a fitted `MarkovChain` as argument `x`.

- k:

  A single positive integer giving the number of trajectories.

- m:

  An integer of at least two giving the number of states observed per
  trajectory (`m - 1` transitions).

- ...:

  Additional arguments passed to `fun`.

- n_workers:

  A single positive integer giving the maximum number of local PSOCK
  workers. Set to one for sequential execution.

## Value

The results of `fun` combined with
[`c()`](https://rdrr.io/r/base/c.html) in simulation order. For example,
if each call returns a numeric vector of length `p`, the combined result
has length `p * k`, with one consecutive block per trajectory. For one
trajectory, the result of `fun` is returned directly.

## Details

Each fit estimates the transition matrix with [`fit()`](fit.md) and
retains the initial distribution from `x`. An error is raised if any
trajectory omits a state, because the fit would have a smaller state
space. Choose a trajectory length that gives adequate coverage of the
states.

Use [`set.seed()`](https://rdrr.io/r/base/Random.html) before calling
this function for reproducible simulations with a fixed worker count.
Workers use the calling session's library paths. For parallel execution,
`fun` and its arguments must be usable in a fresh R session with the
package installed.

## See also

[`sim()`](sim.md), [`fit()`](fit.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
x <- MarkovChain(c(1, 0), P)
set.seed(1)
fit_parallel(x, mttf, k = 4, m = 200, up = 1, n_workers = 1)
#> [1] 5.000000 5.619048 4.346154 4.360000
```
