# Asymptotic Variance of the Mean-Time-to-Recovery Estimator

Computes the asymptotic variance of the plug-in estimator of the mean
time to recovery (MTTR) for a finite-state Markov chain.

## Usage

``` r
var_mttr(x, down)
```

## Arguments

- x:

  An object of class `MarkovChain`.

- down:

  Integer or character vector giving the down (non-operational) states.

## Value

A numeric scalar giving the expected number of discrete-time transitions
before the chain first enters a state outside `down`.

## Details

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

The MTTR \\\psi\\ is returned by [`mttr()`](mttr.md). The asymptotic
variance \\\sigma\_\psi^2\\ is returned by this function.

## See also

[`mttf()`](mttf.md), `var_mttr()`

## Examples

``` r
# Choose the transition probability matrix and the initial law of the MC
a <- c(0.5, 0.5, 0)
P <- matrix(c(0.3, 0.7, 0,
              0.6, 0.1, 0.3,
              0.5, 0,   0.5), nrow = 3, ncol = 3, byrow = TRUE)

# Create a Markov Chain
x <- MarkovChain(a = a, P = P)

# Find Stationary distribution
pi_vec <- stat(x)

# Simulate a single trajectory of length m
m <- 1e4
traj <- sim(m, x)

# Estimate the transition probability matrix
y <- fit(traj)

# Make a list of Markov Chains
L <- MarkovChainList(list(x, y))
L[[1]]
#> Markov Chain
#>   Number of states: 3 
#> 
#>   Initial distribution (a):
#>       1   2   3 
#>     0.5 0.5 0.0 
#> 
#>   Transition matrix (P):
#>         1   2   3
#>     1 0.3 0.7 0.0
#>     2 0.6 0.1 0.3
#>     3 0.5 0.0 0.5
```
