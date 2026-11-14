# BMP Failure Rate of a Markov Chain

Compute the discrete-time BMP failure rate: the probability of first
failure at a given time conditional on survival to the preceding time.

## Usage

``` r
bmp(x, time, up)
```

## Arguments

- x:

  An object of class `MarkovChain`.

- time:

  A nonempty nonnegative integer vector of discrete time points. Values
  must be less than `.Machine$integer.max`; missing and infinite values
  are not allowed.

- up:

  Integer-valued numeric vector of state positions or character vector
  of state names giving the operational states. Missing values are not
  allowed. Duplicates are ignored; empty and full state sets are valid.
  States are used in their order in `x`.

## Value

A numeric vector of the same length and order as `time` containing the
BMP failure rates, including the paper's zero-survival convention.

## Details

Equation (2.7) of Sadek and Limnios (2002) defines
\$\$\lambda(t)=1-R(t)/R(t-1), \qquad t\geq1,\$\$ when \\R(t-1)\>0\\, and
sets \\\lambda(t)=0\\ otherwise. At time zero, \\\lambda(0)=1-R(0)\\ is
the initial probability of being down. Here \\R(t)=a_U
P\_{UU}^t\mathbf{1}\\.

The implementation propagates the conditional distribution among
surviving operational states, normalizing at each step. Its weighted
exit probability equals the reliability ratio above, without subtracting
nearly equal probabilities or forming powers that may underflow. The
initial distribution is treated as fixed. Empty operational sets give
one at time zero and zero thereafter; full operational sets give zero at
every time.

## References

Sadek, A., and Limnios, N. (2002). *Asymptotic properties for maximum
likelihood estimators for reliability and failure rates of Markov
chains*. Communications in Statistics - Theory and Methods, 31(10),
1837–1861.
[doi:10.1081/STA-120014916](https://doi.org/10.1081/STA-120014916) .

## See also

[`rel()`](rel.md), [`ftp()`](ftp.md), [`vcov_bmp()`](vcov_bmp.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
x <- MarkovChain(c(0.4, 0.6), P)
bmp(x, time = 0:3, up = 1L)
#> [1] 0.6 0.2 0.2 0.2
vcov_bmp(x, time = 3L, up = 1L)
#> [1] 0.2666667
```
