# Create a VARBootstrapResult

Create a VARBootstrapResult

## Usage

``` r
create_var_bootstrap_result(
  adjacency_matrices,
  total_effects,
  lags,
  resampled_indices = NULL,
  causal_orders = NULL
)
```

## Arguments

- adjacency_matrices:

  list (length n_sampling); each element is a joined adjacency matrix
  (n_features x n_features\*(1 + lags))

- total_effects:

  array (n_sampling x n_features x n_features\*(1 + lags))

- lags:

  lag order used

- resampled_indices:

  list of residual-index vectors (NULL allowed)

- causal_orders:

  matrix (n_sampling x n_features) (NULL allowed)

## Value

a VARBootstrapResult (list with class attribute)
