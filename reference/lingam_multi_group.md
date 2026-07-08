# Multi-Group Direct LiNGAM

Jointly estimates a Direct LiNGAM model from multiple datasets
("groups") that are assumed to share a common causal order but may have
different structural coefficients (Shimizu 2012).

## Usage

``` r
lingam_multi_group(
  X_list,
  prior_knowledge = NULL,
  apply_prior_knowledge_softly = FALSE,
  reg_method = "adaptive_lasso",
  lambda = "BIC",
  init_method = "ols"
)
```

## Arguments

- X_list:

  A list of numeric matrices or data frames (length \>= 2), one per
  group. Each element must have `n_d` rows (sample size may differ by
  group) and the same number of columns `p` across all groups.

- prior_knowledge:

  Prior knowledge matrix (n_features x n_features) or NULL. Applied
  identically to every group. Same encoding as
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md):
  0: no directed path from x_i to x_j 1: directed path from x_i to x_j
  -1: unknown

- apply_prior_knowledge_softly:

  Whether to apply prior knowledge softly (logical)

- reg_method:

  Regression method for adjacency matrix estimation. "ols": ordinary
  least squares, "lasso": LASSO regression, "adaptive_lasso": adaptive
  LASSO regression (default), "ridge": Ridge regression (robust to
  multicollinearity; does not perform sparse estimation).

- lambda:

  LASSO penalty (lambda) selection. Same choices as
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md).

- init_method:

  Method for estimating the initial weights of adaptive LASSO regression
  ("ols" (default) or "ridge").

## Value

A `MultiGroupLingamResult` object (list) containing:

- `adjacency_matrices`: a named list of adjacency matrices, one per
  group (name = group name). Each follows the same convention as
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md):
  `B[i, j]` is the causal coefficient from variable j to variable i (j
  -\> i).

- `causal_order`: the causal order shared by all groups (integer vector
  of 1-based indices).

## Details

Unlike
[`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md),
this function has no `measure` argument: the multi-group causal-order
search only supports the pwling (pairwise likelihood / entropy
approximation) objective, matching the upstream Python
`MultiGroupDirectLiNGAM`, which does not offer a kernel-based
multi-group search.

For downstream analysis of a single group (total causal effects,
independence tests of residuals, plotting), extract that group as a
plain `LingamResult` with
[`get_group_result()`](https://morimotoosamu.github.io/lingamr/reference/get_group_result.md)
and pass it to the existing single-group functions
([`estimate_total_effect()`](https://morimotoosamu.github.io/lingamr/reference/estimate_total_effect.md),
[`estimate_all_total_effects()`](https://morimotoosamu.github.io/lingamr/reference/estimate_all_total_effects.md),
[`get_error_independence_p_values()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values.md),
[`plot_adjacency()`](https://morimotoosamu.github.io/lingamr/reference/plot_adjacency.md),
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html),
[`tidy()`](https://generics.r-lib.org/reference/tidy.html)); no
multi-group-specific wrappers are provided for these.

## References

S. Shimizu. Joint estimation of linear non-Gaussian acyclic models.
Neurocomputing, 81: 104-107, 2012.

## Examples

``` r
mg <- generate_multi_group_sample()
res <- lingam_multi_group(mg$data_list, reg_method = "ols")
print(res)
#> Multi-Group Direct LiNGAM Result
#>   Groups      : 2 (group1, group2)
#>   Variables   : 6
#>   Causal order (common): x3 -> x0 -> x5 -> x2 -> x4 -> x1
#> 
#> [group1] Adjacency matrix (row = to, col = from):
#>        x0 x1     x2     x3     x4    x5
#> x0  0.000  0  0.000  3.033  0.000 0.000
#> x1  3.237  0  1.965  0.014 -0.034 0.006
#> x2 -0.236  0  0.000  6.112  0.000 0.049
#> x3  0.000  0  0.000  0.000  0.000 0.000
#> x4  7.921  0 -1.063  0.399  0.000 0.018
#> x5  4.016  0  0.000 -0.003  0.000 0.000
#> 
#> [group2] Adjacency matrix (row = to, col = from):
#>       x0 x1     x2     x3    x4     x5
#> x0 0.000  0  0.000  3.504 0.000  0.000
#> x1 2.732  0  2.568  0.083 0.034  0.093
#> x2 0.154  0  0.000  6.322 0.000 -0.024
#> x3 0.000  0  0.000  0.000 0.000  0.000
#> x4 8.483  0 -1.487 -0.110 0.000  0.006
#> x5 4.515  0  0.000 -0.045 0.000  0.000

# Analyze group 1 with the existing single-group tooling
g1 <- get_group_result(res, 1)
estimate_all_total_effects(mg$data_list[[1]], g1, method = "ols")
#>             x0 x1        x2        x3          x4          x5
#> x0  0.00000000  0  0.000000  3.033460  0.00000000  0.00000000
#> x1  2.90911952  0  2.001580 21.058733 -0.03397056  0.10299386
#> x2 -0.03933572  0  0.000000  5.992677  0.00000000  0.04894766
#> x3  0.00000000  0  0.000000  0.000000  0.00000000  0.00000000
#> x4  8.03407606  0 -1.062516 18.276121  0.00000000 -0.03416285
#> x5  4.01586857  0  0.000000 12.179395  0.00000000  0.00000000
```
