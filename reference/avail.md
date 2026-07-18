# Availability of a Markov Chain

Compute the discrete-time availability function of a system whose
degradation/repair dynamics are modeled by a finite-state Markov chain.
Availability at time \\t\\ is defined as \$\$ A(t) = \Pr\\ X_t \in
\mathcal{U} \\ \$\$ where \\X_t\\ is the state of the chain at time
\\t\\ and \\\mathcal{U}\\ is the set of *up* (operational) states.

## Usage

``` r
avail(x, time, up)
```

## Arguments

- x:

  An object of class `MarkovChain`.

- time:

  A positive integer vector giving the discrete time points \\t\\ at
  which availability \\A(t)\\ is evaluated.

- up:

  Integer or character vector giving the up (operational) states.

## Value

A numeric vector of the same length as `time`, where the \\j\\-th
element is the probability that the system is available (in an up state)
at `time[j]`.

## Details

For each time \\t\\, the function computes \$\$ A(t) = a P^t
\mathbf{1}\_{\mathcal{U}}, \$\$ where \\a\\ is the initial distribution,
\\P\\ is the transition matrix, and \\\mathbf{1}\_{\mathcal{U}}\\ is the
indicator vector of the operational states.

Duplicate values in `up` are ignored.

## See also

[`MarkovChain()`](MarkovChain.md), [`Ppower()`](Ppower.md)

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
