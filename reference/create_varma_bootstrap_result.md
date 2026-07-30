# Create a VARMABootstrapResult

Create a VARMABootstrapResult

## Usage

``` r
create_varma_bootstrap_result(
  adjacency_matrices,
  total_effects,
  order,
  resampled_indices = NULL,
  causal_orders = NULL
)
```

## Arguments

- adjacency_matrices:

  list (length n_sampling); each element is a joined matrix
  `cbind(psi_0..psi_p, omega_1..omega_q)` (n_features x n_features\*(1 +
  p + q))

- total_effects:

  array (n_sampling x n_features x n_features\*(1 + p))

- order:

  VARMA order c(p, q) used

- resampled_indices:

  list of residual-index vectors (NULL allowed)

- causal_orders:

  matrix (n_sampling x n_features) (NULL allowed)

## Value

a VARMABootstrapResult (list with class attribute)
