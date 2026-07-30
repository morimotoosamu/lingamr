# Convert a CAMUVResult to a tidy data.frame

Converts the estimated adjacency matrix into a long-format data.frame
with one edge per row, like
[`tidy.LingamResult()`](https://morimotoosamu.github.io/lingamr/reference/tidy.LingamResult.md).
Because CAM-UV is a nonlinear method, `estimate` is 1 for identified
edges (an edge indicator, not a coefficient), and `NA` for variable
pairs suspected to be connected through an unobserved variable
(UCP/UBP); the `NA` rows appear in both directions. Drop them with e.g.
`subset(tidy(x), !is.na(estimate))` if not needed.

## Usage

``` r
# S3 method for class 'CAMUVResult'
tidy(x, threshold = 0, ...)
```

## Arguments

- x:

  The return value of
  [`lingam_camuv()`](https://morimotoosamu.github.io/lingamr/reference/lingam_camuv.md)
  (a `CAMUVResult` object)

- threshold:

  Kept for interface consistency with the other
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) methods
  (default: 0); with 0/1 entries any value in `[0, 1)` returns all
  edges. `NA` entries are always kept.

- ...:

  Unused

## Value

data.frame(from, to, estimate)

## Examples

``` r
# \donttest{
if (requireNamespace("mgcv", quietly = TRUE)) {
  confounded <- generate_camuv_sample(n = 200, seed = 1)
  model <- lingam_camuv(confounded$data)
  tidy(model)
}
#>   from to estimate
#> 1   x0 x1        1
#> 2   x0 x3        1
#> 3   x2 x4        1
#> 4   x2 x5       NA
#> 5   x3 x4       NA
#> 6   x4 x3       NA
#> 7   x5 x2       NA
# }
```
