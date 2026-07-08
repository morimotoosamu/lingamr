# Bottom-Up ParceLiNGAM

A causal ordering method robust against latent confounders. Unlike
[`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md),
which always returns a full causal order, this algorithm searches from
the sink (most downstream) side and stops as soon as an independence
test is rejected. Variables it could not order are returned together as
a single "unresolved block" (suspected to share a latent confounder),
and the corresponding entries of the adjacency matrix are set to `NA`
rather than estimated.

## Usage

``` r
lingam_parce(
  X,
  alpha = 0.1,
  prior_knowledge = NULL,
  independence = "hsic",
  ind_corr = 0.5,
  reg_method = "adaptive_lasso",
  lambda = "BIC",
  init_method = "ols"
)
```

## Arguments

- X:

  Numeric matrix (n_samples x n_features), data frame or matrix

- alpha:

  Significance level for the independence test. `alpha = 0` disables
  rejection entirely, so a full causal order is always returned
  (equivalent in spirit to
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md),
  but using the ParceLiNGAM sink-search direction and regression). Must
  be non-negative.

- prior_knowledge:

  Prior knowledge matrix (n_features x n_features) or NULL. 0: no
  directed path from x_i to x_j 1: directed path from x_i to x_j -1:
  unknown

- independence:

  Independence measure used for the ordering search: "hsic" (default)
  uses the HSIC gamma-approximation test combined across explanatory
  variables via Fisher's method; "fcorr" uses the F-correlation (kernel
  canonical correlation) and rejects based on `ind_corr` instead of a
  p-value.

- ind_corr:

  Threshold on the F-correlation value used only when
  `independence = "fcorr"`: a candidate is rejected once the largest
  F-correlation against any explanatory variable is at or above this
  value. Must be non-negative. Ignored when `independence = "hsic"`.

- reg_method:

  Regression method for adjacency matrix estimation. "ols": ordinary
  least squares, "lasso": LASSO regression, "adaptive_lasso": adaptive
  LASSO regression (default, matches the upstream Python
  implementation's `predict_adaptive_lasso`), "ridge": Ridge regression.

- lambda:

  LASSO penalty (lambda) selection. Same options as
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md):
  "lambda.min", "lambda.1se", "AIC", "BIC" (default), "oracle" (adaptive
  LASSO only).

- init_method:

  Method for estimating the initial weights of adaptive LASSO regression
  ("ols" (default) or "ridge").

## Value

A `ParceLingamResult` object (list) containing:

- `adjacency_matrix`: adjacency matrix B (n_features x n_features).
  **Convention: `B[i, j]` is the causal coefficient from variable j to
  variable i (j -\> i)**, same as
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md).
  Entries between two variables that ended up in the same unresolved
  block are `NA`.

- `causal_order`: a list of integer vectors. Elements of length 1 are
  variables with a fully resolved position; an element of length \> 1
  (at most one, always first) is the unresolved block. Earlier elements
  are more upstream.

- `p_values`: independence-test p-values (or F-correlation values, for
  `independence = "fcorr"`) for each step that successfully placed a
  variable, in the order variables were placed (diagnostic only).

- `independence`: the independence measure used.

## Details

Because HSIC forms full n x n Gram matrices, it is O(n^2) per test;
avoid very large `n` (beyond a few thousand) with
`independence = "hsic"`.

`independence = "fcorr"` rejects based on the raw F-correlation value
(`ind_corr`), not a p-value, so it is not directly comparable to
`alpha`.

[`get_error_independence_p_values_parce()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values_parce.md)
uses the HSIC test rather than the correlation-based test used by
[`get_error_independence_p_values()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values.md)
for `LingamResult` objects.

[`lingam_parce_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce_bootstrap.md)
treats `NA` (unresolved) edges as absent when aggregating, and does not
support
[`get_causal_order_stability()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_order_stability.md)
(see its documentation for details). This function does not expose a
`regressor` or `bw_method` argument, unlike the upstream Python
implementation.

## References

Tashiro, T., Shimizu, S., Hyvarinen, A., and Washio, T. (2014).
ParceLiNGAM: a causal ordering method robust against latent confounders.
Neural Computation, 26(1), 57-83.

## Examples

``` r
confounded <- generate_parce_sample(n = 500, seed = 1)

result <- lingam_parce(confounded$data, reg_method = "ols")
print(result)
#> Bottom-Up ParceLiNGAM Result
#>   Variables : 6
#>   Independence measure: hsic
#>   Causal order: (x2, x3) -> x0 -> x4 -> x5 -> x1
#>   (NA entries in the adjacency matrix = unresolved order / suspected latent confounding)
#> 
#> Adjacency matrix (row = to, col = from):
#>       x0 x1     x2     x3    x4     x5
#> x0 0.000  0 -0.010  0.516 0.000  0.000
#> x1 0.479  0  0.447  0.060 0.025 -0.049
#> x2 0.000  0  0.000     NA 0.000  0.000
#> x3 0.000  0     NA  0.000 0.000  0.000
#> x4 0.497  0 -0.490 -0.001 0.000  0.000
#> x5 0.436  0  0.068  0.023 0.050  0.000

# The variable pair sharing the latent confounder is left unresolved (NA)
result$adjacency_matrix[confounded$confounded_pair, confounded$confounded_pair]
#>    x2 x3
#> x2  0 NA
#> x3 NA  0

# Total effect estimation warns and returns NA for confounded variables
estimate_total_effect_parce(confounded$data, result,
  from_index = confounded$confounded_pair[1], to_index = 1
)
#> Warning: x2 is part of an unresolved causal order (suspected latent confounding); total effect cannot be estimated.
#> [1] NA
```
