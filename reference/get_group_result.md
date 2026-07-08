# Extract a single group's result from a MultiGroupLingamResult

Returns the adjacency matrix and (shared) causal order of one group as a
plain `LingamResult`, so that the existing single-group functions
([`estimate_total_effect()`](https://morimotoosamu.github.io/lingamr/reference/estimate_total_effect.md),
[`estimate_all_total_effects()`](https://morimotoosamu.github.io/lingamr/reference/estimate_all_total_effects.md),
[`get_error_independence_p_values()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values.md),
[`plot_adjacency()`](https://morimotoosamu.github.io/lingamr/reference/plot_adjacency.md),
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html),
[`tidy()`](https://generics.r-lib.org/reference/tidy.html)) can be
applied to it directly.

## Usage

``` r
get_group_result(x, group)
```

## Arguments

- x:

  A `MultiGroupLingamResult`, as returned by
  [`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md).

- group:

  Group name (character) or 1-based group index (integer).

## Value

A `LingamResult` object (list) with `adjacency_matrix` and
`causal_order`, identical in shape to the return value of
[`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md).

## Examples

``` r
mg <- generate_multi_group_sample()
res <- lingam_multi_group(mg$data_list, reg_method = "ols")
g1 <- get_group_result(res, 1)
class(g1)
#> [1] "LingamResult"
```
