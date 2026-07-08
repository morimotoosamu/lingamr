# Causal-order search for HighDimDirectLiNGAM

Causal-order search for HighDimDirectLiNGAM

## Usage

``` r
high_dim_causal_order(X, J, K, alpha)
```

## Arguments

- X:

  data matrix (n_samples x n_features)

- J:

  assumed largest in-degree

- K:

  moment degree

- alpha:

  pruning cutoff coefficient

## Value

integer vector, 1-based causal order (upstream-most first)
