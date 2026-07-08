# Collapse an ImputationBootstrapResult into a BootstrapResult

[`bootstrap_with_imputation()`](https://morimotoosamu.github.io/lingamr/reference/bootstrap_with_imputation.md)'s
result carries an extra `n_repeats` dimension (one causal-structure
estimate per imputed dataset per bootstrap iteration), which the
existing bootstrap analysis functions
([`get_probabilities()`](https://morimotoosamu.github.io/lingamr/reference/get_probabilities.md),
[`get_causal_direction_counts()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_direction_counts.md),
[`get_directed_acyclic_graph_counts()`](https://morimotoosamu.github.io/lingamr/reference/get_directed_acyclic_graph_counts.md),
[`get_causal_order_stability()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_order_stability.md),
[`tidy()`](https://generics.r-lib.org/reference/tidy.html)) do not
expect. This collapses that dimension by aggregating the `n_repeats`
adjacency matrices of each iteration into one, producing a regular
[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md)-style
`BootstrapResult`.

## Usage

``` r
as_bootstrap_result(x, aggregate = c("median", "mean"))
```

## Arguments

- x:

  An `ImputationBootstrapResult`, as returned by
  [`bootstrap_with_imputation()`](https://morimotoosamu.github.io/lingamr/reference/bootstrap_with_imputation.md).

- aggregate:

  How to aggregate across the `n_repeats` dimension: `"median"`
  (default) or `"mean"`.

## Value

A `BootstrapResult` (see
[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md))
with `total_effects = NULL` (total effects are not computed by
[`bootstrap_with_imputation()`](https://morimotoosamu.github.io/lingamr/reference/bootstrap_with_imputation.md));
calling
[`get_total_causal_effects()`](https://morimotoosamu.github.io/lingamr/reference/get_total_causal_effects.md)
on it raises the usual "not computed" error.

## Examples

``` r
set.seed(1)
sample6 <- generate_lingam_sample_6(n = 300, seed = 1)
X <- sample6$data
X$x5[sample.int(nrow(X), size = 30)] <- NA

if (requireNamespace("mice", quietly = TRUE)) {
  res <- bootstrap_with_imputation(X,
    n_sampling = 5L, n_repeats = 3L, seed = 42, verbose = FALSE
  )
  bs <- as_bootstrap_result(res, aggregate = "median")
  get_probabilities(bs)

  # get_total_causal_effects() is not available: total effects were never computed
  tryCatch(get_total_causal_effects(bs), error = function(e) conditionMessage(e))
}
#> [1] "This BootstrapResult has no total effects. It was created with compute_total_effects = FALSE; re-run lingam_direct_bootstrap() with compute_total_effects = TRUE to use get_total_causal_effects()."
```
