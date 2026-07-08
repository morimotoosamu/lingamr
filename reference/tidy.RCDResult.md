# Convert an RCDResult to a tidy data.frame

Converts the estimated adjacency matrix into a long-format data.frame
with one edge per row, like
[`tidy.LingamResult()`](https://morimotoosamu.github.io/lingamr/reference/tidy.LingamResult.md).
`NA` entries of the adjacency matrix (variable pairs suspected to share
a latent confounder; marked in both directions) are kept as rows with
`estimate = NA` so they remain visible; drop them with e.g.
`subset(tidy(x), !is.na(estimate))` if not needed.

## Usage

``` r
# S3 method for class 'RCDResult'
tidy(x, threshold = 0, ...)
```

## Arguments

- x:

  The return value of
  [`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)
  (an `RCDResult` object)

- threshold:

  Coefficients with an absolute value at or below this are not treated
  as edges (default: 0). `NA` entries are always kept.

- ...:

  Unused

## Value

data.frame(from, to, estimate)

## Examples

``` r
confounded <- generate_rcd_sample(n = 300, seed = 1)
model <- lingam_rcd(confounded$data)
tidy(model)
#>   from to  estimate
#> 1   x0 x2 0.8095226
#> 2   x0 x4 1.0153709
#> 3   x1 x0 1.1158054
#> 4   x2 x4        NA
#> 5   x3 x0 0.9887483
#> 6   x4 x2        NA
#> 7   x5 x1 0.5880495
#> 8   x5 x3 0.4486123
```
