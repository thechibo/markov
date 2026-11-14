# Simulate Markov Chain Trajectories

Generate one or more independent sample paths from a finite-state,
time-homogeneous Markov chain.

## Usage

``` r
sim(m, x, n = 1)
```

## Arguments

- m:

  A single positive integer giving the number of observations per
  trajectory, including the initial state.

- x:

  A `MarkovChain` object.

- n:

  A single positive integer giving the number of trajectories.

## Value

If `n = 1`, an ordered factor of length `m` with one level per state in
chain order. If `n > 1`, an `n` by `m` character matrix of state labels,
with one trajectory per row and observations in time order across
columns. Pass a matrix with at least two columns directly to
[`fit()`](fit.md) to pool its within-trajectory transitions.

## Details

Each trajectory starts with a fresh draw from the initial distribution
`x$a`. Each subsequent state is drawn from the row of `x$P` for the
current state. Simulating multiple paths with a fixed random seed gives
the same draws as successive single-path calls with that seed.

## See also

[`MarkovChain()`](MarkovChain.md), [`fit()`](fit.md),
[`stat()`](stat.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
x <- MarkovChain(a = c(0.5, 0.5), P = P)
set.seed(1)
sim(m = 5, x = x)
#> [1] 2 2 2 1 1
#> Levels: 1 < 2
paths <- sim(m = 20, x = x, n = 10)
fit(paths)
#> Markov Chain
#>   Number of states: 2 
#> 
#>   Initial distribution (a):
#>       1   2 
#>     0.2 0.8 
#> 
#>   Transition matrix (P):
#>               1         2
#>     1 0.7634409 0.2365591
#>     2 0.2783505 0.7216495
```
