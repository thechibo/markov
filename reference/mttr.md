# Mean Time to Recovery

Computes the mean time to recovery of a discrete-time, finite-state
Markov chain. Recovery is defined as the first time the chain leaves the
set of non-operational states specified by `down`.

## Usage

``` r
mttr(x, down)
```

## Arguments

- x:

  An object of class `MarkovChain`.

- down:

  Nonempty integer-valued numeric vector of state positions, or
  character vector of state names, giving the non-operational states.
  Entries must be distinct and select a proper subset of the states.
  Input order is ignored; states are used in their order in `x`.

## Value

A numeric scalar giving the expected number of discrete-time transitions
before the chain first enters a state outside `down`.

## Details

These functions use [`mttf()`](mttf.md) and [`var_mttf()`](var_mttf.md)
with `down` as the selected state set. The non-operational subchain must
be transient, so that \\I_d-P\_{\mathcal{D},\mathcal{D}}\\ is
invertible. The asymptotic variance additionally requires an irreducible
chain and treats the initial distribution as fixed.

The initial probabilities on `down` are not renormalized. Paths starting
outside `down` have recovery time zero; to compute the mean conditional
on starting down, supply an initial distribution supported on `down`.

Let \\\mathcal{D}\\ denote the set of non-operational states,
\\P\_{\mathcal{D},\mathcal{D}}\\ be the submatrix of transition
probabilities between non-operational states, \\d=\|\mathcal{D}\|\\, and
\\N=(I_d-P\_{\mathcal{D},\mathcal{D}})^{-1}\\. The mean time to recovery
is computed as \$\$ \psi = a\_{\mathcal{D}}^\mathsf{T} N \mathbf{1}\_d,
\$\$ where \\a\_{\mathcal{D}}\\ is the restriction of the initial
distribution to the non-operational states and \\\mathbf{1}\_d\\ is a
vector of ones.

Let \\\Sigma\_{\mathcal{D}}\\ be the submatrix of \\\Sigma\_{\mathbf
p}\\ corresponding to transitions between non-operational states. Then:
\$\$ \sqrt{n} \left\\ \widehat{\psi}\_n-\psi \right\\
\mathrel{\mathop{\longrightarrow}\limits^{d}} \mathcal{N} \left( 0,
\sigma\_\psi^2 \right), \$\$ where \$\$ \sigma\_\psi^2 = D\_\psi
\Sigma\_{\mathcal{D}} D\_\psi^{\mathsf T}, \$\$ and \$\$ D\_\psi =
\left( a\_{\mathcal{D}}^{\mathsf T}N \right) \mathbin{\otimes} \left(
\mathbf{1}\_d^{\mathsf T}N^{\mathsf T} \right). \$\$

The MTTR \\\psi\\ is returned by this function. The asymptotic variance
\\\sigma\_\psi^2\\ is returned by [`var_mttr()`](var_mttr.md).

## References

Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026). *Matrix
asymptotic calculus for plug-in maximum likelihood estimators in finite
Markov chains*. arXiv preprint arXiv:2607.17161.
<https://arxiv.org/abs/2607.17161>. See Proposition 11 for the MTTR
differential and asymptotic variance.

Sadek, A., and Limnios, N. (2002). *Asymptotic properties for maximum
likelihood estimators for reliability and failure rates of Markov
chains*. Communications in Statistics - Theory and Methods, 31(10),
1837–1861.
[doi:10.1081/STA-120014916](https://doi.org/10.1081/STA-120014916) .

## See also

[`mttf()`](mttf.md), [`var_mttr()`](var_mttr.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
dimnames(P) <- list(c("working", "failed"), c("working", "failed"))
x <- MarkovChain(c(working = 0, failed = 1), P)
mttr(x, down = "failed")
#> [1] 3.333333
var_mttr(x, down = "failed")
#> [1] 64.81481
```
