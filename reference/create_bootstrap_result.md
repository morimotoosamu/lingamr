# Create a BootstrapResult

Create a BootstrapResult

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

  matrix (n_sampling x n_features). Each row is the causal order of one
  sample.

## Value

BootstrapResult (list with class attribute)
