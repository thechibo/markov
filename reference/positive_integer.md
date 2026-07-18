# Validate or Coerce Positive Integer Values

Utilities for validating and coercing positive integer-valued numeric
vectors.

`is_positive_integer()` tests whether every element of `m` is finite,
integer-valued, greater than or equal to one, and representable as a
base R integer.

`as_positive_integer()` performs the same validation and converts `m` to
an integer vector. It raises an error when validation fails.

## Usage

``` r
is_positive_integer(m)

as_positive_integer(m)
```

## Arguments

- m:

  An object to validate or coerce.

## Value

`is_positive_integer()` returns a single logical value.

`as_positive_integer()` returns `m` as an integer vector, or raises an
error if `m` is not a valid positive integer vector.

## Details

A value is considered a positive integer when it:

- is numeric;

- has positive length;

- contains no missing or infinite values;

- contains only whole numbers;

- is greater than or equal to one; and

- does not exceed `.Machine$integer.max`.

Zero, negative numbers, fractional values, `NA`, `NaN`, `Inf`, and
non-numeric objects are rejected.

## Examples

``` r
is_positive_integer(3)
#> [1] TRUE
is_positive_integer(c(1, 2, 3))
#> [1] TRUE
is_positive_integer(c(1, 2.5))
#> [1] FALSE
is_positive_integer(0)
#> [1] FALSE

as_positive_integer(3)
#> [1] 3
as_positive_integer(c(1, 2, 3))
#> [1] 1 2 3
```
