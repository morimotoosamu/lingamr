# Get a one-row summary of a ResitResult

Like
[`glance.LingamResult()`](https://morimotoosamu.github.io/lingamr/reference/glance.LingamResult.md),
with an additional `regressor` column giving the regressor label used
for the internal nonlinear regressions.

## Usage

``` r
# S3 method for class 'ResitResult'
glance(x, ...)
```

## Arguments

- x:

  The return value of
  [`lingam_resit()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit.md)
  (a `ResitResult` object)

- ...:

  Unused

## Value

A one-row data.frame(n_variables, n_edges, regressor, causal_order)

## Examples

``` r
# \donttest{
if (requireNamespace("mgcv", quietly = TRUE)) {
  nonlinear <- generate_resit_sample(n = 300, seed = 1)
  model <- lingam_resit(nonlinear$data)
  glance(model)
}
#>   n_variables n_edges regressor         causal_order
#> 1           4       4       gam x0 -> x1 -> x2 -> x3
# }
```
