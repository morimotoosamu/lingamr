# Median-heuristic kernel width for HSIC

Uses only the first 100 rows (not a random subsample) to keep the O(n^2)
pairwise-distance computation cheap, matching the upstream
implementation exactly. Multivariate input is measured by the squared
Euclidean distance between rows, as in the upstream
`get_kernel_width()`.

## Usage

``` r
hsic_kernel_width(x)
```

## Arguments

- x:

  numeric vector or matrix (n x d)

## Value

kernel width (scalar)
