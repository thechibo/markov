# Mean Time to Failure

Mean Time to Failure

## Usage

``` r
mttf(x, up)

var_mttf(x, up)

var_mttf_sl(x, up)

deriv_psi(x, up)
```

## Arguments

- x:

  an object of class `MarkovChain`.

- up:

  integer or character. The up states.

- ...:

  extra arguments.

## Value

numeric. The Mean Time to Failure.

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
