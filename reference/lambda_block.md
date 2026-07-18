# Row-Wise Asymptotic Covariance Block of Transition Probabilities

Computes the \\m \times m\\ covariance block \$\$ \Lambda_i =
\operatorname{diag}(p\_{i\bullet}) - p\_{i\bullet}p\_{i\bullet}^{\mathsf
T}, \$\$ associated with row \\i\\ of a finite-state Markov chain
transition matrix \\P\\. Here, \\p\_{i\bullet} = (p\_{i1}, \ldots,
p\_{im})^{\mathsf T}\\ is the vector of transition probabilities from
state \\i\\.

## Usage

``` r
lambda_block(p)
```

## Arguments

- p:

  Numeric probability row vector \\p\_{i\bullet}\\. Should be
  non-negative and sum to 1 (row of a row-stochastic transition matrix).

## Value

A numeric matrix of dimension \\s \times s\\, where \\s =
\mathrm{length}(p)\\ is the number of states.

## Details

This matrix is the covariance matrix of a single multinomial trial with
category-probability vector \\p\_{i\bullet}\\. It also forms the
\\i\\-th diagonal block of the asymptotic covariance matrix of the
empirical transition-probability MLE \\\widehat{\mathbf P}(n)\\.

The \\(j,k)\\ entry of the matrix is \$\$ \Lambda_i(j,k) =
p\_{ij}\left(\delta\_{jk} - p\_{ik}\right), \qquad j,k = 1,\ldots,m,
\$\$ where \\\delta\_{jk}\\ is the Kronecker delta.

## See also

[`vcov_p()`](vcov_p.md)
