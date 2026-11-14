# Mean Time to Failure

Computes the mean time to failure of a discrete-time, finite-state
Markov chain. Failure is defined as the first time the chain leaves the
set of operational states specified by `up`.

## Usage

``` r
mttf(x, up)
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

A numeric scalar giving the expected number of discrete-time transitions
before the chain first enters a state outside `up`.

## Details

The operational subchain must be transient, so that
\\I_r-P\_{\mathcal{U},\mathcal{U}}\\ is invertible. The asymptotic
variance additionally requires an irreducible chain and treats the
initial distribution as fixed.

Let \\\mathcal{U}\\ denote the set of operational states,
\\P\_{\mathcal{U},\mathcal{U}}\\ be the submatrix of transition
probabilities between operational states, \\r=\|\mathcal{U}\|\\, and
\\N=(I_r-P\_{\mathcal{U},\mathcal{U}})^{-1}\\. The mean time to failure
is computed as \$\$ \psi = a\_{\mathcal{U}}^\mathsf{T} N \mathbf{1}\_r,
\$\$ where \\a\_{\mathcal{U}}\\ is the restriction of the initial
distribution to the operational states and \\\mathbf{1}\_r\\ is a vector
of ones.

Let \\\Sigma\_{\mathcal{U}}\\ be the submatrix of \\\Sigma\_{\mathbf
p}\\ corresponding to transitions between operational states. Then: \$\$
\sqrt{n} \left\\ \widehat{\psi}\_n-\psi \right\\
\mathrel{\mathop{\longrightarrow}\limits^{d}} \mathcal{N} \left( 0,
\sigma\_\psi^2 \right), \$\$ where \$\$ \sigma\_\psi^2 = D\_\psi
\Sigma\_{\mathcal{U}} D\_\psi^{\mathsf T}, \$\$ and \$\$ D\_\psi =
\left( a\_{\mathcal{U}}^{\mathsf T}N \right) \mathbin{\otimes} \left(
\mathbf{1}\_r^{\mathsf T}N^{\mathsf T} \right). \$\$

The MTTF \\\psi\\ is returned by this function. The asymptotic variance
\\\sigma\_\psi^2\\ is returned by [`var_mttf()`](var_mttf.md).

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

[`avail()`](avail.md) for the time-dependent availability function.

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
x <- MarkovChain(c(1, 0), P)
mttf(x, up = 1)
#> [1] 5
var_mttf(x, up = 1)
#> [1] 166.6667
```
