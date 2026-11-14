# Test a Specified Markov Transition Probability

Test the null hypothesis that a single transition probability equals
`p0`, against a two-sided alternative, using the asymptotic test of
Anderson and Goodman (1957).

## Usage

``` r
test_prob(x, ij, p0)
```

## Arguments

- x:

  An atomic vector or factor of states in time order, or an atomic
  matrix with one independent trajectory per row and common observation
  times in columns. Each path must have at least two observations;
  matrices must have at least one row. Missing values are not allowed.
  Validated by the internal `validate_trajectories()` helper.

- ij:

  An atomic vector or factor of length two, with no dimensions or
  missing values. Its first entry is the source state and its second is
  the destination state. Both must occur in `x`. Entries are state
  values, not positions in the state space; equal entries test a
  self-transition.

- p0:

  A single finite numeric probability strictly between zero and one.
  Boundary nulls are excluded because the asymptotic null variance
  vanishes there.

## Value

An object of class `htest`, printed using the standard R test format,
with components `statistic` (named `X-squared`), `parameter` (one degree
of freedom), `p.value`, `estimate`, `null.value`, `alternative`
(`"two.sided"`), `method`, and `data.name`. Additional components
`observed` and `expected` contain the two-category observed and null
expected transition counts.

## Details

Section 3.1 (pp. 96–97) of Anderson and Goodman (1957) derives tests for
specified transition probabilities from their multinomial asymptotic
law. For a single entry, aggregate destinations into state \\j\\ and its
complement. Writing \\N_i\\ for the number of departures from state
\\i\\, \\N\_{ij}\\ for transitions from \\i\\ to \\j\\, and \\\widehat
p\_{ij}=N\_{ij}/N_i\\, the two-category version of Equation (3.1) is
\$\$X^2 = \frac{N_i(\widehat p\_{ij}-p_0)^2}{p_0(1-p_0)}.\$\$ Under the
null, this has an asymptotic chi-squared distribution with one degree of
freedom. No continuity correction is applied. The final observation is
not counted as a departure; an error is raised if the source state has
no observed departures.

The model is a time-homogeneous, first-order finite Markov chain. For a
single long trajectory, the approximation requires adequate visits to
the source state, as under the ergodic setting discussed in Section 5.
It is not an exact binomial test conditional on the random departure
count. For a matrix, rows are independent paths with a common transition
matrix; counts are pooled across paths, excluding each path's final
observation from departures. The Section 3.1 approximation applies as
the number of independent paths increases with their common length
fixed. A warning is issued if either null expected count, \\N_i p_0\\ or
\\N_i(1-p_0)\\, is below five. The result reports a p-value rather than
imposing a significance level; reject at level \\\alpha\\ when the
p-value is below \\\alpha\\.

## References

Anderson, T. W., and Goodman, L. A. (1957). Statistical Inference about
Markov Chains. *The Annals of Mathematical Statistics*, 28(1), 89–110.
[doi:10.1214/aoms/1177707039](https://doi.org/10.1214/aoms/1177707039) .
See Section 3.1 for specified-probability tests and Section 5 for a
single long chain.

## See also

[`fit()`](fit.md),
[`stats::chisq.test()`](https://rdrr.io/r/stats/chisq.test.html),
[`stats::prop.test()`](https://rdrr.io/r/stats/prop.test.html)

## Examples

``` r
x <- rep(c("a", "a", "b", "a", "b", "b"), 20)
test_prob(x, ij = c("a", "b"), p0 = 0.5)
#> 
#>  Anderson-Goodman test of a transition probability
#> 
#> data:  x
#> X-squared = 6.6667, df = 1, p-value = 0.009823
#> alternative hypothesis: true p[a, b] is not equal to 0.5
#> sample estimates:
#>   p[a, b] 
#> 0.6666667 
#> 
```
