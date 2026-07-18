# Kronecker Delta Vector

Creates a vector of length `n` with a one at position `i` and zeros
elsewhere.

## Usage

``` r
cron(i, n)
```

## Arguments

- i:

  Position of the nonzero element.

- n:

  Length of the vector.

## Value

A numeric vector of length `n`.

## Examples

``` r
cron(2, 4)
#> [1] 0 1 0 0
```
