# BootstrapResult を作成

BootstrapResult を作成

## Usage

``` r
create_bootstrap_result(
  adjacency_matrices,
  total_effects,
  resampled_indices = NULL,
  causal_orders = NULL
)
```

## Arguments

- adjacency_matrices:

  array (n_sampling x n_features x n_features)

- total_effects:

  array (n_sampling x n_features x n_features)

- resampled_indices:

  list of index vectors

- causal_orders:

  matrix (n_sampling x n_features)。各行が1標本の因果順序。

## Value

BootstrapResult (list with class attribute)
