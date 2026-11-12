# Asymptotic Variance–Covariance of the Transition Probability Vector

Constructs the block-diagonal asymptotic variance–covariance matrix for
the row-wise vectorization of an empirical transition matrix from an
irreducible finite-state Markov chain.

## Usage

``` r
vcov_p(x, ...)

# S3 method for class 'MarkovChain'
vcov_p(x, ...)

# S3 method for class 'StochMat'
vcov_p(x, ...)

# Default S3 method
vcov_p(x, ...)
```

## Arguments

- x:

  A `MarkovChain` or `StochMat` object representing an \\s \times s\\
  row-stochastic transition matrix. The Markov chain must be
  irreducible.

- ...:

  Additional arguments passed to methods.

## Value

A numeric block-diagonal matrix of dimension \\s^2 \times s^2\\, where
\\s\\ is the number of states.

## Details

Let \\\widehat{\mathbf P}(n)\\ be the empirical transition matrix
estimated from a trajectory of length \\n\\, and let \\\widehat{\mathbf
p}(n)\\ be its row-wise vectorization. Under irreducibility, the
empirical transition matrix is strongly consistent and asymptotically
normal: \$\$ \sqrt{n} \left\\ \widehat{\mathbf p}(n) - \mathbf p
\right\\ \xrightarrow\[n\to\infty\]{d} \mathbf W\_{\mathbf p}, \qquad
\mathbf W\_{\mathbf p} \sim \mathcal N\_{s^2} \left( \mathbf 0\_{s^2},
\Sigma\_{\mathbf p} \right). \$\$

This function calculates the matrix \\\Sigma\_{\mathbf p}\\.
Specifically, if \\\boldsymbol{\pi} = (\pi_1,\ldots,\pi_s)^{\mathsf T}\\
is the unique stationary distribution, then the asymptotic covariance
matrix is \$\$ \Sigma\_{\mathbf p} = \begin{pmatrix} \pi_1^{-1}\Lambda_1
& 0 & \cdots & 0 \\ 0 & \pi_2^{-1}\Lambda_2 & \cdots & 0 \\ \vdots &
\vdots & \ddots & \vdots \\ 0 & 0 & \cdots & \pi_s^{-1}\Lambda_s
\end{pmatrix}, \$\$

where \$\$ \Lambda_i = \operatorname{diag}(p\_{i\bullet}) -
p\_{i\bullet}p\_{i\bullet}^{\mathsf T}, \qquad i = 1,\ldots,s. \$\$

## References

Gavrilopoulos, G., Trevezas, S., and Votsi, I. (2026). *Matrix
asymptotic calculus for plug-in maximum likelihood estimators in finite
Markov chains*. arXiv preprint arXiv:2607.17161.
<https://arxiv.org/abs/2607.17161>. See Theorem 1 for the row-wise
transition covariance.

Sadek, A., and Limnios, N. (2002). *Asymptotic properties for maximum
likelihood estimators for reliability and failure rates of Markov
chains*. Communications in Statistics - Theory and Methods, 31(10),
1837–1861.
[doi:10.1081/STA-120014916](https://doi.org/10.1081/STA-120014916) . See
Proposition 1 (p. 1846) for the asymptotic covariance matrix of the
transition probabilities.

## Examples

``` r
P <- StochMat(matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE))
Sigma <- vcov_p(P)
dim(Sigma) # Four transition probabilities in row-wise order
#> [1] 4 4
```
