# Kernel-based mutual information

Dispatches to the incomplete-Cholesky low-rank path for n above the
low-rank threshold (matching the kappa/sigma switch in
[`search_causal_order_kernel()`](https://morimotoosamu.github.io/lingamr/reference/search_causal_order_kernel.md));
below the threshold it calls the exact path unchanged.

## Usage

``` r
mutual_information_kernel(x1, x2, param)
```

## Arguments

- x1:

  Variable 1

- x2:

  Variable 2

- param:

  Parameter vector (kappa, sigma)

## Value

Mutual information
