# Get a one-row summary of a CAMUVResult

Like
[`glance.LingamResult()`](https://morimotoosamu.github.io/lingamr/reference/glance.LingamResult.md),
but without a causal order (CAM-UV does not estimate one). `n_edges`
counts non-`NA` edges only, and `n_confounded_pairs` counts the variable
pairs whose adjacency-matrix entries are `NA` (suspected unobserved
causal/backdoor path).

## Usage

``` r
# S3 method for class 'CAMUVResult'
glance(x, ...)
```

## Arguments

- x:

  The return value of
  [`lingam_camuv()`](https://morimotoosamu.github.io/lingamr/reference/lingam_camuv.md)
  (a `CAMUVResult` object)

- ...:

  Unused

## Value

A one-row data.frame(n_variables, n_edges, n_confounded_pairs,
regressor)

## Examples

``` r
# \donttest{
if (requireNamespace("mgcv", quietly = TRUE)) {
  confounded <- generate_camuv_sample(n = 200, seed = 1)
  model <- lingam_camuv(confounded$data)
  glance(model)
}
#>   n_variables n_edges n_confounded_pairs regressor
#> 1           6       3                  2       gam
# }
```
