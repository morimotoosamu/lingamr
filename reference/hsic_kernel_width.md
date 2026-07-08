# Median-heuristic kernel width for HSIC

Uses only the first 100 points (not a random subsample) to keep the
O(n^2) pairwise-distance computation cheap, matching the upstream
implementation exactly.

## Usage

``` r
hsic_kernel_width(x)
```

## Arguments

- x:

  numeric vector

## Value

kernel width (scalar)
