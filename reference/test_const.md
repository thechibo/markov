# Test Constancy of Markov Transition Probabilities across Time

Test whether independent Markov trajectories share a transition matrix
that is constant across observation times, using the Pearson test in
Section 3.2 of Anderson and Goodman (1957), Equation (3.8).

## Usage

``` r
test_const(x)
```

## Arguments

- x:

  An atomic matrix with one independent trajectory per row and common
  observation times in columns. At least two rows and three columns are
  required. Missing values are not allowed. Validated by the internal
  `validate_trajectories()` helper. A single vector is not sufficient
  for this test.

## Value

An object of class `htest`, with components `statistic` (named
`X-squared`), `parameter` (degrees of freedom), `p.value`, `method`, and
`data.name`. Additional components are `estimate` (the pooled transition
matrix), and `observed` and `expected` (arrays indexed by source state,
destination state, and transition time). State order is first appearance
when reading paths row by row.

## Details

Write \\T = \mathrm{ncol}(x)-1\\ for the number of transition times and
\\m\\ for the number of observed states. Let \\N\_{ij}(t)\\ count paths
making transition \\i\to j\\ at time \\t\\, and let \\N_i(t-1)=\sum_j
N\_{ij}(t)\\. Under the null, the pooled estimate is \$\$\widehat
p\_{ij} = \frac{\sum_t N\_{ij}(t)}{\sum_t N_i(t-1)}.\$\$ The statistic
is \$\$X^2 = \sum\_{i,j,t} \frac{(N\_{ij}(t)-N_i(t-1)\widehat
p\_{ij})^2} {N_i(t-1)\widehat p\_{ij}}.\$\$ Its asymptotic null
distribution is chi-squared with \\m(m-1)(T-1)\\ degrees of freedom. The
alternative allows a different transition matrix at each time. Paths
must be independent, but observations within a path need not be
independent. The approximation is for increasing numbers of paths with
fixed \\m\\ and \\T\\, positive transition probabilities, and growing
origin counts at every transition time.

This implementation requires at least two observed states, a positive
origin count at every transition time for every state, and positive
pooled counts for every ordered state pair. Empty margins cause an error
rather than an automatic change to the paper's degrees of freedom. A
warning is issued when any expected cell count is below five. No
continuity correction is applied. Reject at a chosen significance level
when the p-value is below that level. This tests time homogeneity, not
equality of the initial and stationary distributions.

## References

Anderson, T. W., and Goodman, L. A. (1957). Statistical Inference about
Markov Chains. *The Annals of Mathematical Statistics*, 28(1), 89–110.
[doi:10.1214/aoms/1177707039](https://doi.org/10.1214/aoms/1177707039) .
Section 3.2, especially Equations (3.6) and (3.8).

## See also

[`test_prob()`](test_prob.md), [`fit()`](fit.md),
[`stats::chisq.test()`](https://rdrr.io/r/stats/chisq.test.html)

## Examples

``` r
# Repeated copies of all binary paths of length three.
paths <- as.matrix(expand.grid(start = c("a", "b"),
                              middle = c("a", "b"), end = c("a", "b")))
x <- paths[rep(seq_len(nrow(paths)), each = 10), ]
test_const(x)
#> 
#>  Anderson-Goodman test of constant transition probabilities
#> 
#> data:  x
#> X-squared = 0, df = 2, p-value = 1
#> sample estimates:
#>     to
#> from   a   b
#>    a 0.5 0.5
#>    b 0.5 0.5
#> 
```
