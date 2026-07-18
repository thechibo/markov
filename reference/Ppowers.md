# Successive Powers of a Transition Matrix

Computes successive powers of a transition matrix by repeated matrix
multiplication.

## Usage

``` r
Ppowers(P, k)
```

## Arguments

- P:

  A square numeric transition matrix.

- k:

  A positive integer giving the highest matrix power to compute.

## Value

A list of `k + 1` matrices. Element `j + 1` contains \\P^j\\, beginning
with the identity matrix \\P^0\\.

## Details

Unlike [`Ppower()`](Ppower.md), this function does not use an
eigendecomposition. Instead, each power is obtained by multiplying the
preceding power by `P`:

\$\$P^j = P^{j-1}P.\$\$

Repeated multiplication can be useful when all powers from \\P^0\\
through \\P^k\\ are required.

The current implementation expects `k` to be at least one.

## See also

[`Ppower()`](Ppower.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)

powers <- Ppowers(P, k = 3)

powers[[1]] # P^0
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
powers[[2]] # P^1
#>      [,1] [,2]
#> [1,]  0.8  0.2
#> [2,]  0.3  0.7
powers[[4]] # P^3
#>       [,1]  [,2]
#> [1,] 0.650 0.350
#> [2,] 0.525 0.475
```
