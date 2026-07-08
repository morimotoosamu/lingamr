# Get a one-row summary of a ParceLingamResult

Like
[`glance.LingamResult()`](https://morimotoosamu.github.io/lingamr/reference/glance.LingamResult.md).
`n_edges` counts non-`NA` edges only, and `n_na_entries` counts the
adjacency-matrix entries left `NA` (unresolved order / suspected latent
confounding). Unresolved blocks in the causal order are shown in
parentheses, as in the print method.

## Usage

``` r
# S3 method for class 'ParceLingamResult'
glance(x, ...)
```

## Arguments

- x:

  The return value of
  [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)
  (a `ParceLingamResult` object)

- ...:

  Unused

## Value

A one-row data.frame(n_variables, n_edges, n_na_entries, causal_order)

## Examples

``` r
dat <- generate_parce_sample(n = 500, seed = 42)
model <- lingam_parce(dat$data)
glance(model)
#>   n_variables n_edges n_na_entries                     causal_order
#> 1           6       7            0 x2 -> x3 -> x0 -> x5 -> x1 -> x4
```
