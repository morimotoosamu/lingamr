# Get a one-row summary of an RCDResult

Like
[`glance.LingamResult()`](https://morimotoosamu.github.io/lingamr/reference/glance.LingamResult.md),
but without a causal order (RCD does not estimate one). `n_edges` counts
non-`NA` edges only, and `n_confounded_pairs` counts the variable pairs
whose adjacency-matrix entries are `NA` (suspected shared latent
confounder).

## Usage

``` r
# S3 method for class 'RCDResult'
glance(x, ...)
```

## Arguments

- x:

  The return value of
  [`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)
  (an `RCDResult` object)

- ...:

  Unused

## Value

A one-row data.frame(n_variables, n_edges, n_confounded_pairs)

## Examples

``` r
confounded <- generate_rcd_sample(n = 300, seed = 1)
model <- lingam_rcd(confounded$data)
glance(model)
#>   n_variables n_edges n_confounded_pairs
#> 1           6       6                  1
```
