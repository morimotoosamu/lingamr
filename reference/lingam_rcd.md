# RCD (Repetitive Causal Discovery)

A causal discovery method robust against latent confounders. Unlike
[`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
or
[`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md),
RCD does not attempt to recover a full or partial causal order. Instead,
it repeatedly extracts each variable's **ancestor set** by scanning
variable subsets of increasing size
([`extract_ancestors()`](https://morimotoosamu.github.io/lingamr/reference/extract_ancestors.md)),
narrows each ancestor set down to direct **parents**
([`extract_parents()`](https://morimotoosamu.github.io/lingamr/reference/extract_parents.md)),
and finally tests remaining parent-free pairs for a shared latent
confounder
([`extract_vars_sharing_confounders()`](https://morimotoosamu.github.io/lingamr/reference/extract_vars_sharing_confounders.md)).
Pairs found to share a latent confounder are marked `NA` in the
adjacency matrix rather than estimated.

## Usage

``` r
lingam_rcd(
  X,
  max_explanatory_num = 2L,
  cor_alpha = 0.01,
  ind_alpha = 0.01,
  shapiro_alpha = 0.01,
  MLHSICR = FALSE,
  independence = "hsic",
  ind_corr = 0.5
)
```

## Arguments

- X:

  Numeric matrix (n_samples x n_features), data frame or matrix

- max_explanatory_num:

  Maximum number of explanatory variables considered when searching for
  ancestors (i.e. the search scans variable subsets of size up to
  `max_explanatory_num + 1`). Larger values increase statistical power
  but grow combinatorially in cost. Must be an integer of 1 or more.

- cor_alpha:

  Significance level for the Pearson correlation tests used throughout
  the algorithm (ancestor-subset screening, parent extraction,
  confounder-pair detection). Must be non-negative.

- ind_alpha:

  Significance level for the HSIC independence test (used when
  `independence = "hsic"`). Must be non-negative.

- shapiro_alpha:

  Significance level for the Shapiro-Wilk non-Gaussianity test used when
  screening candidate ancestor subsets. Must be non-negative.

- MLHSICR:

  If `TRUE`, falls back to HSIC-sum-minimizing regression (instead of
  OLS) when the OLS residual is not independent of the explanatory
  variables and more than one explanatory variable is present.
  Substantially increases computation time.

- independence:

  Independence measure used for the sink search: "hsic" (default) uses
  the HSIC gamma-approximation test; "fcorr" uses the F-correlation
  (kernel canonical correlation) and rejects based on `ind_corr` instead
  of a p-value.

- ind_corr:

  Threshold on the F-correlation value, used only when
  `independence = "fcorr"`. Must be non-negative. Ignored when
  `independence = "hsic"`.

## Value

An `RCDResult` object (list) containing:

- `adjacency_matrix`: adjacency matrix B (n_features x n_features).
  **Convention: `B[i, j]` is the causal coefficient from variable j to
  variable i (j -\> i)**, same as
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md).
  Entries between two variables found to share a latent confounder are
  `NA`.

- `ancestors_list`: a list of length `n_features`; element `i` is the
  sorted integer vector of variables found to be ancestors of variable
  `i` (possibly empty). Unlike
  [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md),
  there is no `causal_order`: RCD estimates ancestor relations directly
  rather than a total or partial order.

## Details

The algorithm has three stages: (1)
[`extract_ancestors()`](https://morimotoosamu.github.io/lingamr/reference/extract_ancestors.md)
grows each variable's ancestor set by repeatedly scanning variable
subsets; (2)
[`extract_parents()`](https://morimotoosamu.github.io/lingamr/reference/extract_parents.md)
narrows ancestor sets down to direct parents; (3)
[`extract_vars_sharing_confounders()`](https://morimotoosamu.github.io/lingamr/reference/extract_vars_sharing_confounders.md)
tests remaining parent-free pairs for a shared latent confounder. `NA`
entries in `adjacency_matrix` mean the corresponding pair is suspected
to share a latent confounder, not that no relationship was estimated.

`max_explanatory_num` controls both statistical power and computational
cost: stage 1 scans `choose(n_features, k)` subsets for each subset size
`k` up to `max_explanatory_num + 1`, and each subset requires several
HSIC tests (each O(n^2) in the sample size when
`independence = "hsic"`). Cost grows quickly with both the number of
variables and `n`.

`MLHSICR = TRUE` replaces the OLS residual in the independence check
with a residual obtained by directly minimizing the sum of HSIC
statistics between the residual and each explanatory variable via
numerical optimization (`stats::optim(method = "L-BFGS-B")`). This can
recover independence in cases where OLS cannot, but requires
re-optimizing for every candidate subset where the OLS residual fails,
and is therefore substantially slower.

The Shapiro-Wilk test
([`stats::shapiro.test()`](https://rdrr.io/r/stats/shapiro.test.html))
used for the non-Gaussianity check is limited to `n` between 3 and 5000.
For `n` above 5000, a deterministic evenly-spaced subsample of 5000
observations is tested instead (same policy as
[`test_residual_normality()`](https://morimotoosamu.github.io/lingamr/reference/test_residual_normality.md)),
so results remain reproducible without touching the RNG state. This
subsampling has no effect when `n <= 5000`.

This function does not expose a `bw_method` argument (kernel widths are
always the median heuristic; see
[`hsic_kernel_width()`](https://morimotoosamu.github.io/lingamr/reference/hsic_kernel_width.md)),
unlike some upstream implementations.
[`lingam_rcd_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd_bootstrap.md)
does not support
[`get_causal_order_stability()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_order_stability.md),
since RCD has no causal order.

## References

Maeda, T. N. and Shimizu, S. (2020). RCD: Repetitive causal discovery of
linear non-Gaussian acyclic models with latent confounders. AISTATS
2020, PMLR 108: 735-745.

## Examples

``` r
confounded <- generate_rcd_sample(n = 300, seed = 1)

result <- lingam_rcd(confounded$data)
print(result)
#> RCD Result
#>   Variables : 6
#> 
#> Ancestor sets:
#>   M(x0) = {x1, x3, x5}
#>   M(x1) = {x5}
#>   M(x2) = {x0, x1, x3, x5}
#>   M(x3) = {x5}
#>   M(x4) = {x0, x1, x3, x5}
#>   M(x5) = {}
#> 
#>   (NA entries in the adjacency matrix = suspected shared latent confounder)
#> 
#> Adjacency matrix (row = to, col = from):
#>       x0    x1 x2    x3 x4    x5
#> x0 0.000 1.116  0 0.989  0 0.000
#> x1 0.000 0.000  0 0.000  0 0.588
#> x2 0.810 0.000  0 0.000 NA 0.000
#> x3 0.000 0.000  0 0.000  0 0.449
#> x4 1.015 0.000 NA 0.000  0 0.000
#> x5 0.000 0.000  0 0.000  0 0.000

# The variable pair sharing the latent confounder is left NA
result$adjacency_matrix[confounded$confounded_pair, confounded$confounded_pair]
#>    x2 x4
#> x2  0 NA
#> x4 NA  0

# Total effect estimation warns and returns NA for confounded variables
estimate_total_effect_rcd(confounded$data, result,
  from_index = confounded$confounded_pair[1], to_index = 1
)
#> Warning: x0 is an ancestor of x2 according to ancestors_list; the requested direction (x2 -> x0) is inconsistent with the estimated ancestor relations.
#> Warning: x2 is part of a suspected latent confounder pair; total effect cannot be estimated.
#> [1] NA
```
