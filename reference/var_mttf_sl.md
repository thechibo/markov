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

  Nonempty integer-valued numeric vector of state positions, or
  character vector of state names, giving the operational states.
  Entries must be distinct and select a proper subset of the states.
  Input order is ignored; states are used in their order in `x`.

## Value

`var_mttf_sl()` returns a numeric scalar containing the asymptotic
variance. `deriv_psi()` returns a numeric vector of length
`length(up)^2`, containing the derivatives with respect to the
operational transition block in row-wise order, with operational states
ordered as in `x`.

## Details

Unlike [`var_mttf()`](var_mttf.md), which uses a matrix and
Kronecker-product implementation, this function follows the
corresponding scalar or component-wise expressions. It is used for
internal testing only.

## References

Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026). *Matrix
asymptotic calculus for plug-in maximum likelihood estimators in finite
Markov chains*. arXiv preprint arXiv:2607.17161.
<https://arxiv.org/abs/2607.17161>. See Proposition 11 for the MTTF
differential and asymptotic variance.

Sadek, A., and Limnios, N. (2002). *Asymptotic properties for maximum
likelihood estimators for reliability and failure rates of Markov
chains*. Communications in Statistics - Theory and Methods, 31(10),
1837–1861.
[doi:10.1081/STA-120014916](https://doi.org/10.1081/STA-120014916) .

## See also

[`mttf()`](mttf.md), [`var_mttf()`](var_mttf.md)
