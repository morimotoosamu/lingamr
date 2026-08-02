# Precompute the per-variable parts of the HSIC gamma test

Computes everything
[`hsic_test_gamma()`](https://morimotoosamu.github.io/lingamr/reference/hsic_test_gamma.md)
derives from one argument alone (kernel width, centered Gram matrix, and
the mean term `mu`), so callers that test many pairs sharing a variable
can reuse the O(n^2) Gram computation. Validation (n \>= 6, no NA,
constant input) matches
[`hsic_test_gamma()`](https://morimotoosamu.github.io/lingamr/reference/hsic_test_gamma.md)
exactly.

## Usage

``` r
hsic_precompute(x)
```

## Arguments

- x:

  numeric vector or matrix (n x d)

## Value

list(n, is_const, Kc = centered Gram matrix, mu = mean term); `Kc`/`mu`
are NULL when `is_const` is TRUE
