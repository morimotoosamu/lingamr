# Bootstrap for RESIT

Bootstrap for RESIT

## Usage

``` r
lingam_resit_bootstrap(
  X,
  n_sampling,
  regressor = "gam",
  alpha = 0.01,
  prior_knowledge = NULL,
  seed = NULL,
  verbose = TRUE,
  parallel = FALSE,
  n_cores = NULL
)
```

## Arguments

- X:

  Numeric matrix (n_samples x n_features)

- n_sampling:

  Number of bootstrap iterations

- regressor:

  Nonlinear regressor, passed to
  [`lingam_resit()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit.md)

- alpha:

  Significance level of the HSIC pruning test, passed to
  [`lingam_resit()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit.md)

- prior_knowledge:

  Prior-knowledge matrix, passed to
  [`lingam_resit()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit.md)

- seed:

  Random seed (NULL allowed)

- verbose:

  Whether to display progress (logical)

- parallel:

  Whether to use parallel processing (logical)

- n_cores:

  Number of cores to use (integer, NULL allowed)

## Value

A `BootstrapResult` (list); see
[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md)
for the query helpers that operate on it
([`get_probabilities()`](https://morimotoosamu.github.io/lingamr/reference/get_probabilities.md),
[`get_causal_direction_counts()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_direction_counts.md),
[`get_directed_acyclic_graph_counts()`](https://morimotoosamu.github.io/lingamr/reference/get_directed_acyclic_graph_counts.md),
[`get_causal_order_stability()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_order_stability.md)).

## Details

**`total_effects` is always `NULL`**: RESIT is a nonlinear method for
which total causal effects are undefined, so there is no
`compute_total_effects` argument and
[`get_total_causal_effects()`](https://morimotoosamu.github.io/lingamr/reference/get_total_causal_effects.md)
raises its usual "no total effects" error. (The Python implementation
instead stores an all-zero total-effects array; storing nothing is
deliberate, so the zeros cannot be mistaken for estimated effects.)

**`causal_orders` is populated** (RESIT estimates a full causal order),
so
[`get_causal_order_stability()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_order_stability.md)
works with the returned object.

Each iteration runs `O(ncol(X)^2)` nonlinear regressions plus HSIC
tests, so a bootstrap is substantially slower than the linear variants;
keep `n_sampling` modest. When `parallel = TRUE`, a user-supplied
`regressor` function is serialized to the PSOCK workers; any package it
uses must be referenced with the `pkg::fun` form inside the function
body.

## Examples

``` r
# \donttest{
if (requireNamespace("mgcv", quietly = TRUE)) {
  nonlinear <- generate_resit_sample(n = 300, seed = 1)

  bs <- lingam_resit_bootstrap(nonlinear$data,
    n_sampling = 3L,
    seed = 42
  )
  get_probabilities(bs)
}
#> Bootstrap: 3 iterations, RESIT (sequential)
#>   iteration 1 / 3
#> Completed in 1.1 seconds.
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    0    0    0
#> [2,]    1    0    0    0
#> [3,]    1    1    0    0
#> [4,]    0    0    1    0
# }
```
