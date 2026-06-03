# BootstrapResult を作成

BootstrapResult を作成

## Usage

``` r
create_bootstrap_result(
  adjacency_matrices,
  total_effects,
  resampled_indices = NULL
)
```

## Arguments

- adjacency_matrices:

  array (n_sampling x n_features x n_features)

- total_effects:

  array (n_sampling x n_features x n_features)

- resampled_indices:

  list of index vectors

## Value

BootstrapResult (list with class attribute)
