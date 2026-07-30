# Normalize HSIC input to an (n, d) numeric matrix

Vectors become one-column matrices so that all downstream computations
can work row-wise, matching the upstream `X.reshape(-1, 1)` behavior.

## Usage

``` r
as_hsic_matrix(x)
```

## Arguments

- x:

  numeric vector or matrix

## Value

numeric matrix (n x d)
