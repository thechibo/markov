# Availability of a Markov chain

Compute the discrete-time availability function of a system whose
degradation/repair dynamics are modeled by a finite-state Markov chain.
Availability at time \\t\\ is defined as \$\$ A(t) = \Pr\\ X_t \in
\mathcal{U} \\ \$\$ where \\X_t\\ is the state of the chain at time
\\t\\ and \\\mathcal{U}\\ is the set of *up* (operational) states.

Equivalently, for each time \\t\\ in `time`, the function computes \$\$
A(t) = a P^t \mathbf{1}\_{\mathcal{U}}, \$\$ where \\P = x\$P\\ is the
transition matrix and \\\mathbf{1}\_{\mathcal{U}}\\ is the indicator
vector of up states.

## Usage

``` r
avail(x, time, up)
```

## Arguments

- x:

  An object of class `MarkovChain`.

- time:

  Integer vector. Discrete time steps \\t\\ at which to evaluate the
  availability \\A(t)\\.

## Value

A numeric vector of the same length as `time`, where the \\j\\-th
element is the probability that the system is available (in an up state)
at time `time[j]`.

## See also

[`MarkovChain()`](MarkovChain.md), [`mat_gamma()`](mat_gamma.md)

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
pi_vec <- mc_stat(x)

# Simulate a single trajectory of length m
m <- 1e4
traj <- mc_sim(m, x)

# Estimate the transition probability matrix
y <- mc_estim(traj)

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
