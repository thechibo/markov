# Reference Implementation of the Asymptotic Variance of MTTF

Computes the asymptotic variance of the mean-time-to-failure estimator
using a component-wise derivative calculation.

## Usage

``` r
var_mttf_sl(x, up)

deriv_psi(x, up)
```

## Arguments

- x:

  An object of class `MarkovChain`.

- up:

  Integer or character vector giving the up (operational) states.

## Value

A numeric scalar containing the asymptotic variance of the MTTF
estimator.

## Details

Unlike [`var_mttf()`](var_mttf.md), which uses a matrix and
Kronecker-product implementation, this function follows the
corresponding scalar or component-wise expressions. It is used for
internal testing only.

## References

A. Sadek and N. Limnios (2002). *Asymptotic properties of estimators for
nonhomogeneous Markov chains*. Statistics & Probability Letters, 60(4),
1839–1850. See Theorem 8 (p. 1852) for the MTTF.

## See also

[`mttf()`](mttf.md), [`var_mttf()`](var_mttf.md)
