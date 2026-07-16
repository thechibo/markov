# Row-wise asymptotic covariance block of transition probabilities

Compute the \\m \times m\\ block \$\$\Gamma_i =
\mathrm{diag}(p\_{i\bullet}) - p\_{i\bullet} p\_{i\bullet}^\top\$\$
corresponding to a single row \\p\_{i\bullet}\\ of a Markov chain
transition matrix \\P\\, where \\p\_{i\bullet}\\ denotes the row vector
of transition probabilities from state \\i\\.

## Usage

``` r
lambda_block(p)
```

## Arguments

- p:

  Numeric probability row vector \\p\_{i\bullet}\\. Should be
  non-negative and sum to 1 (row of a row-stochastic transition matrix).

## Value

A numeric matrix of dimension \\m \times m\\, where \\m =
\mathrm{length}(p)\\.
