# Bootstrap for RCD

Bootstrap for RCD

## Usage

``` r
lingam_rcd_bootstrap(
  X,
  n_sampling,
  max_explanatory_num = 2L,
  cor_alpha = 0.01,
  ind_alpha = 0.01,
  shapiro_alpha = 0.01,
  MLHSICR = FALSE,
  independence = "hsic",
  ind_corr = 0.5,
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

- max_explanatory_num:

  Maximum number of explanatory variables, passed to
  [`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)

- cor_alpha:

  Significance level for correlation tests, passed to
  [`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)

- ind_alpha:

  Significance level for the HSIC independence test, passed to
  [`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)

- shapiro_alpha:

  Significance level for the non-Gaussianity test, passed to
  [`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)

- MLHSICR:

  Whether to use MLHSICR regression, passed to
  [`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)

- independence:

  Independence measure, passed to
  [`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)

- ind_corr:

  F-correlation rejection threshold, passed to
  [`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)

- seed:

  Random seed (NULL allowed)

- verbose:

  Whether to display progress (logical)

- parallel:

  Whether to use parallel processing (logical)

- n_cores:

  Number of cores to use (integer, NULL allowed)

- compute_total_effects:

  Whether to also estimate total causal effects for every ancestor pair
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

**Total effects are computed only for ancestor pairs**, driven by each
iteration's `ancestors_list` (unlike
[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md)
and
[`lingam_parce_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce_bootstrap.md),
which loop over all variable pairs). For a pair `(from, to)` with `from`
in `to`'s ancestor set, the total effect is obtained via
[`calculate_total_effect()`](https://morimotoosamu.github.io/lingamr/reference/calculate_total_effect.md)
(summing products of adjacency-matrix coefficients along every directed
path). If `from`'s row in the adjacency matrix contains `NA` (it shares
a latent confounder with some other variable), the effect is set to `NA`
for that iteration.

**`NA` (confounded) edges are treated as absent when aggregating.** Both
the adjacency matrix and the total-effect matrix have `NA` replaced by
`0` before being stored in the returned `BootstrapResult`, matching the
policy used by
[`lingam_parce_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce_bootstrap.md).

**`causal_orders` is not populated**: RCD has no causal order (see
[`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)).
As a result,
[`get_causal_order_stability()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_order_stability.md)
cannot be used with a `BootstrapResult` returned by this function.

RCD's `fit` step is HSIC-heavy and can be slow per iteration, especially
with `MLHSICR = TRUE`; keep `n_sampling` modest in examples.

## Examples

``` r
# \donttest{
confounded <- generate_rcd_sample(n = 300, seed = 1)

bs <- lingam_rcd_bootstrap(confounded$data,
  n_sampling = 5L,
  seed = 42
)
#> Bootstrap: 5 iterations, RCD (sequential)
#>   iteration 1 / 5
#> Completed in 3.1 seconds.
get_probabilities(bs)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]  0.0  0.4    0  0.0    0  0.4
#> [2,]  0.0  0.0    0  0.0    0  0.2
#> [3,]  0.2  0.0    0  0.4    0  0.2
#> [4,]  0.0  0.0    0  0.0    0  0.0
#> [5,]  1.0  0.0    0  0.0    0  0.0
#> [6,]  0.0  0.0    0  0.0    0  0.0
# }
```
