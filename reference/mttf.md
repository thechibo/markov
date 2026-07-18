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

  Integer or character vector giving the up (operational) states.

## Value

A numeric scalar giving the expected number of discrete-time transitions
before the chain first enters a state outside `up`.

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

The MTTF \\\psi\\ is returned by this function. The asymptotic variance
\\\sigma\_\psi^2\\ is returned by [`var_mttf()`](var_mttf.md).

## See also

[`avail()`](avail.md) for the time-dependent availability function.

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
