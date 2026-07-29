# Bootstrap for Bottom-Up ParceLiNGAM

Bootstrap for Bottom-Up ParceLiNGAM

## Usage

``` r
lingam_parce_bootstrap(
  X,
  n_sampling,
  prior_knowledge = NULL,
  alpha = 0.1,
  independence = "hsic",
  ind_corr = 0.5,
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

- alpha:

  Significance level, passed to
  [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)

- independence:

  Independence measure, passed to
  [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)

- ind_corr:

  F-correlation rejection threshold, passed to
  [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)

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

  Number of cores to use (integer, NULL allowed)

- compute_total_effects:

  Whether to also estimate total causal effects for every variable pair
  on each bootstrap iteration (logical, default `TRUE`).

## Value

A `BootstrapResult` (list); see
[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md)
for the query helpers that operate on it
([`get_probabilities()`](https://morimotoosamu.github.io/lingamr/reference/get_probabilities.md),
[`get_causal_direction_counts()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_direction_counts.md),
[`get_directed_acyclic_graph_counts()`](https://morimotoosamu.github.io/lingamr/reference/get_directed_acyclic_graph_counts.md),
[`get_total_causal_effects()`](https://morimotoosamu.github.io/lingamr/reference/get_total_causal_effects.md)).

## Details

**Total effects are path sums, not regression estimates.** Each
iteration's total-effect matrix is built from
[`calculate_total_effect()`](https://morimotoosamu.github.io/lingamr/reference/calculate_total_effect.md)
(summing products of adjacency-matrix coefficients along every directed
path), matching the upstream Python implementation's bootstrap method
(`estimate_total_effect2`). If a variable's row in the adjacency matrix
contains `NA` (it is part of an unresolved block), all of its outgoing
total effects are set to `NA` for that iteration, since its causal
parents cannot be identified.

**`NA` (unresolved) edges are treated as absent when aggregating.** Both
the adjacency matrix and the total-effect matrix have `NA` replaced by
`0` before being stored in the returned `BootstrapResult`, matching the
numpy comparison semantics used by the upstream implementation (where
`np.abs(nan) > threshold` evaluates to `FALSE`). This means, for
example,
[`get_probabilities()`](https://morimotoosamu.github.io/lingamr/reference/get_probabilities.md)
reports the confounded pair's edge probability as the fraction of
resamples in which the order happened to resolve, not as `NA`.

**`causal_orders` is not populated** (unlike
[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md)):
ParceLiNGAM's causal order can include an unresolved block, which does
not fit the fixed-length integer-vector format `causal_orders` requires.
As a result,
[`get_causal_order_stability()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_order_stability.md)
cannot be used with a `BootstrapResult` returned by this function.

## Examples

``` r
# \donttest{
confounded <- generate_parce_sample(n = 500, seed = 1)

bs <- lingam_parce_bootstrap(confounded$data,
  n_sampling = 10L,
  reg_method = "ols",
  seed = 42
)
#> Bootstrap: 10 iterations, method=ols (sequential)
#>   iteration 1 / 10
#>   iteration 10 / 10
#> Completed in 15.3 seconds.
get_probabilities(bs)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]  0.0  0.2  0.4  0.4  0.0  0.0
#> [2,]  0.4  0.0  0.4  0.4  0.1  0.3
#> [3,]  0.0  0.0  0.0  0.0  0.0  0.0
#> [4,]  0.0  0.1  0.1  0.0  0.0  0.0
#> [5,]  0.6  0.6  0.6  0.6  0.0  0.4
#> [6,]  0.4  0.3  0.4  0.4  0.2  0.0
# }
```
