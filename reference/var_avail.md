# Variance of availability

Compute the variance of the system availability at discrete time \\n\\
for a finite-state ergodic Markov chain, using the asymptotic covariance
structure described in Theorem 5 of Sadek and Limnios (2002).

## Usage

``` r
var_avail(x, time, up)
```

## Arguments

- x:

  An object of class `MarkovChain`.

- time:

  Integer scalar. Discrete time step \\time\\ at which the variance of
  availability is evaluated.

- up:

  Integer vector giving the indices of up (operational) states.

## Value

A numeric scalar: the variance of the availability \\A(n)\\ at time
\\n\\ as defined in Sadek and Limnios (2002).

## References

A. Sadek and N. Limnios (2002). *Asymptotic properties of estimators for
nonhomogeneous Markov chains*. Statistics & Probability Letters, 60(4),
1839–1850.

## See also

[`avail()`](avail.md), [`mat_gamma()`](mat_gamma.md),
[`mat_sigma()`](mat_sigma.md)
