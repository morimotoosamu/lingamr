# Estimate the adjacency matrix by causal order using the n \<= p route

Estimate the adjacency matrix by causal order using the n \<= p route

## Usage

``` r
estimate_adjacency_matrix_high_dim_np(X, causal_order)
```

## Arguments

- X:

  original-scale data matrix

- causal_order:

  integer vector, 1-based causal order

## Value

adjacency matrix B (n_features x n_features)
