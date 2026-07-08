# Bootstrap for Multi-Group Direct LiNGAM

Bootstrap for Multi-Group Direct LiNGAM

## Usage

``` r
lingam_multi_group_bootstrap(
  X_list,
  n_sampling,
  prior_knowledge = NULL,
  apply_prior_knowledge_softly = FALSE,
  reg_method = "adaptive_lasso",
  lambda = "BIC",
  init_method = "ols",
  seed = NULL,
  verbose = TRUE,
  parallel = FALSE,
  n_cores = NULL,
  compute_total_effects = TRUE
)
```

## Arguments

- X_list:

  A list of numeric matrices or data frames (length \>= 2), one per
  group. Same requirements as
  [`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md).

- n_sampling:

  Number of bootstrap iterations

- prior_knowledge:

  Prior knowledge matrix (NULL allowed). Applied to every group, same as
  [`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md).

- apply_prior_knowledge_softly:

  Apply prior knowledge softly (logical)

- reg_method:

  Regression method ("ols", "lasso", "adaptive_lasso", "ridge")

- lambda:

  Lambda selection ("lambda.min", "lambda.1se", "AIC", "BIC", "oracle")

- init_method:

  Method for estimating the initial weights of adaptive LASSO regression
  ("ols" or "ridge")

- seed:

  Random seed (NULL allowed)

- verbose:

  Whether to display progress (logical)

- parallel:

  Whether to use parallel processing (logical)

- n_cores:

  Number of cores to use (integer, NULL allowed). When `NULL`, the
  number of cores is limited to a maximum of 2 for safety. Ignored when
  `parallel = FALSE`.

- compute_total_effects:

  Whether to also compute total causal effects for every variable pair
  on each bootstrap iteration (logical, default `TRUE`). Set to `FALSE`
  to skip it when only edge/order stability is needed.

## Value

A `MultiGroupBootstrapResult`: a named list (one element per group) of
`BootstrapResult` objects (see
[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md)),
with class `"MultiGroupBootstrapResult"`.

## Details

Each element of the returned list is a regular `BootstrapResult`, so the
existing single-group bootstrap functions
([`get_probabilities()`](https://morimotoosamu.github.io/lingamr/reference/get_probabilities.md),
[`get_causal_direction_counts()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_direction_counts.md),
[`get_directed_acyclic_graph_counts()`](https://morimotoosamu.github.io/lingamr/reference/get_directed_acyclic_graph_counts.md),
[`get_total_causal_effects()`](https://morimotoosamu.github.io/lingamr/reference/get_total_causal_effects.md),
[`get_causal_order_stability()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_order_stability.md),
[`plot_bootstrap_probabilities()`](https://morimotoosamu.github.io/lingamr/reference/plot_bootstrap_probabilities.md),
[`tidy()`](https://generics.r-lib.org/reference/tidy.html)) all work by
extracting a group with `result[[group_name]]` or `result[[i]]`,
mirroring the upstream Python API (which returns a list of
`BootstrapResult` per group).

**Total effects use path products, not regression.** Each bootstrap
iteration's total-effect matrix is the sum of path-coefficient products
over the DAG defined by that iteration's adjacency matrix (matching the
upstream Python MultiGroup bootstrap), which is a different method from
[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md)'s
regression-based
[`estimate_all_total_effects()`](https://morimotoosamu.github.io/lingamr/reference/estimate_all_total_effects.md).

**On iteration failures:** as in
[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md),
each iteration is wrapped in
[`tryCatch()`](https://rdrr.io/r/base/conditions.html); a failing
iteration is skipped with a warning, and only if every iteration fails
is an error raised.

**On reproducibility:** same policy as
[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md).
During parallel execution, L'Ecuyer parallel random number streams via
[`parallel::clusterSetRNGStream()`](https://rdrr.io/r/parallel/RngStream.html)
are used. Results are reproducible given the same `seed` and same
`n_cores`, but they do not numerically match the results of sequential
execution (`parallel = FALSE`). If you need results that exactly match
the sequential version, use `parallel = FALSE`.

## Examples

``` r
mg <- generate_multi_group_sample()

bs <- lingam_multi_group_bootstrap(mg$data_list,
  n_sampling = 10L,
  reg_method = "ols",
  seed = 42
)
#> Multi-group bootstrap: 10 iterations, 2 groups, method=ols (sequential)
#>   iteration 1 / 10
#>   iteration 10 / 10
#> Completed in 0.2 seconds.
get_probabilities(bs[[1]])
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]  0.0  0.0  0.4    1  0.0  0.0
#> [2,]  1.0  0.0  1.0    1  0.6  0.8
#> [3,]  0.6  0.0  0.0    1  0.0  0.4
#> [4,]  0.0  0.0  0.0    0  0.0  0.0
#> [5,]  1.0  0.4  1.0    1  0.0  0.7
#> [6,]  1.0  0.2  0.6    1  0.3  0.0

# \donttest{
bs_par <- lingam_multi_group_bootstrap(mg$data_list,
  n_sampling = 30L,
  reg_method = "ols",
  seed = 42,
  parallel = TRUE,
  n_cores = 2L
)
#> Multi-group bootstrap: 30 iterations, 2 groups, method=ols (parallel, 2 cores)
#> Completed in 0.6 seconds.
# }
```
