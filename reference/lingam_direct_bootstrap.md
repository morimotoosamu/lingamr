# Bootstrap for Direct LiNGAM

Bootstrap for Direct LiNGAM

## Usage

``` r
lingam_direct_bootstrap(
  X,
  n_sampling,
  prior_knowledge = NULL,
  apply_prior_knowledge_softly = FALSE,
  measure = "pwling",
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

- X:

  Numeric matrix (n_samples x n_features)

- n_sampling:

  Number of bootstrap iterations

- prior_knowledge:

  Prior knowledge matrix (NULL allowed)

- apply_prior_knowledge_softly:

  Apply prior knowledge softly (logical)

- measure:

  Independence measure ("pwling" or "kernel")

- reg_method:

  Regression method ("ols", "lasso", "adaptive_lasso", "ridge")

- lambda:

  Lambda selection ("lambda.min", "lambda.1se", "AIC", "BIC","oracle")

- init_method:

  Method for estimating the initial weights of adaptive LASSO regression
  ("ols" or "ridge"). Same as the argument of the same name in
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md).

- seed:

  Random seed (NULL allowed)

- verbose:

  Whether to display progress (logical)

- parallel:

  Whether to use parallel processing (logical). When `TRUE`, each
  bootstrap iteration is distributed across multiple cores.

- n_cores:

  Number of cores to use (integer, NULL allowed). When `NULL`, the
  number of cores is limited to a maximum of 2 for safety. Ignored when
  `parallel = FALSE`.

- compute_total_effects:

  Whether to also estimate total causal effects for every variable pair
  on each bootstrap iteration (logical, default `TRUE`). For the
  lasso-family regression methods this roughly doubles iteration cost
  (an additional regression per downstream variable beyond the
  adjacency-matrix fit). Set to `FALSE` to skip it when only edge/order
  stability is needed
  ([`get_probabilities()`](https://morimotoosamu.github.io/lingamr/reference/get_probabilities.md),
  [`get_causal_direction_counts()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_direction_counts.md),
  [`get_directed_acyclic_graph_counts()`](https://morimotoosamu.github.io/lingamr/reference/get_directed_acyclic_graph_counts.md),
  [`get_causal_order_stability()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_order_stability.md));
  in that case
  [`get_total_causal_effects()`](https://morimotoosamu.github.io/lingamr/reference/get_total_causal_effects.md)
  errors if called on the result.

## Value

BootstrapResult (list)

## Details

When `parallel = TRUE` is specified, iterations are distributed across a
socket cluster created by
[`parallel::makePSOCKcluster()`](https://rdrr.io/r/parallel/makeCluster.html).
The cluster is always released via
[`on.exit()`](https://rdrr.io/r/base/on.exit.html), whether the process
finishes normally or an error occurs.

**On iteration failures:** each bootstrap iteration is run inside a
[`tryCatch()`](https://rdrr.io/r/base/conditions.html). If an iteration
errors (e.g. a resample produces near-singular columns), a warning
identifying the failed iteration is issued and that iteration is
excluded from the result instead of aborting the entire run. The
returned `BootstrapResult` reflects however many iterations actually
succeeded; an error is raised only if every iteration fails.

**On reproducibility:** During parallel execution, L'Ecuyer parallel
random number streams via
[`parallel::clusterSetRNGStream()`](https://rdrr.io/r/parallel/RngStream.html)
are used. Results are reproducible given the same `seed` and same
`n_cores`, but they do not numerically match the results of sequential
execution (`parallel = FALSE`). If you need results that exactly match
the sequential version, use `parallel = FALSE`.

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

# Fast example with OLS
bs <- lingam_direct_bootstrap(LiNGAM_sample_1000$data,
  n_sampling = 10L,
  reg_method = "ols",
  seed = 42
)
#> Bootstrap: 10 iterations, method=ols (sequential)
#>   iteration 1 / 10
#>   iteration 10 / 10
#> Completed in 0.1 seconds.
get_probabilities(bs)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]  0.0  0.1  0.5  0.9  0.1  0.0
#> [2,]  0.9  0.0  0.9  0.9  0.4  0.5
#> [3,]  0.5  0.1  0.0  0.9  0.1  0.4
#> [4,]  0.1  0.1  0.1  0.0  0.1  0.1
#> [5,]  0.9  0.6  0.9  0.9  0.0  0.5
#> [6,]  1.0  0.5  0.6  0.9  0.5  0.0

# \donttest{
# With LASSO (requires glmnet)
bs_lasso <- lingam_direct_bootstrap(LiNGAM_sample_1000$data,
  n_sampling = 30L,
  seed = 42
)
#> Bootstrap: 30 iterations, method=adaptive_lasso (sequential)
#>   iteration 1 / 30
#>   iteration 10 / 30
#>   iteration 20 / 30
#>   iteration 30 / 30
#> Completed in 0.8 seconds.

# Parallel execution on 2 cores
bs_par <- lingam_direct_bootstrap(LiNGAM_sample_1000$data,
  n_sampling = 30L,
  seed = 42,
  parallel = TRUE,
  n_cores = 2L
)
#> Bootstrap: 30 iterations, method=adaptive_lasso (parallel, 2 cores)
#> Completed in 1.4 seconds.
# }
```
