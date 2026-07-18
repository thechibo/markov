# Asymptotic Variance of the Mean-Time-to-Failure Estimator

Computes the asymptotic variance of the plug-in estimator of the mean
time to failure (MTTF) for a finite-state Markov chain.

## Usage

``` r
var_mttf(x, up)
```

## Arguments

- x:

  An object of class `MarkovChain`.

- up:

  Integer or character vector giving the up (operational) states.

## Value

A numeric scalar containing the asymptotic variance \\\sigma\_\psi^2\\
of the MTTF estimator.

## Details

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

The MTTF \\\psi\\ is returned by [`mttf()`](mttf.md). The asymptotic
variance \\\sigma\_\psi^2\\ is returned by this function.

## References

GSI

A. Sadek and N. Limnios (2002). *Asymptotic properties of estimators for
nonhomogeneous Markov chains*. Statistics & Probability Letters, 60(4),
1839–1850. See Theorem 8 (p. 1852) for the MTTF.

## See also

[`mttf()`](mttf.md)
