# Variance of availability

Compute the variance of the system availability at discrete time \\n\\
for a finite-state ergodic Markov chain.

## Usage

``` r
vcov_avail(x, time, up)
```

## Arguments

- x:

  An object of class `MarkovChain`.

- time:

  Integer scalar. Discrete time step \\time\\ at which the variance of
  availability is evaluated.

- up:

  Integer or character vector giving the up (operational) states.

## Value

A numeric scalar: the variance of the availability \\A(n)\\ at time
\\n\\ as defined in Sadek and Limnios (2002).

## References

A. Sadek and N. Limnios (2002). *Asymptotic properties of estimators for
nonhomogeneous Markov chains*. Statistics & Probability Letters, 60(4),
1839–1850. See Theorem 5.

## See also

[`avail()`](avail.md), [`vcov_p()`](vcov_p.md),
[`vcov_pk()`](vcov_pk.md)
