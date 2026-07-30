# Convert a ResitResult to a tidy data.frame

Converts the estimated 0/1 adjacency matrix into a long-format
data.frame with one edge per row, like
[`tidy.LingamResult()`](https://morimotoosamu.github.io/lingamr/reference/tidy.LingamResult.md).
Because RESIT is a nonlinear method, `estimate` is always 1 (an edge
indicator, not a coefficient).

## Usage

``` r
# S3 method for class 'ResitResult'
tidy(x, threshold = 0, ...)
```

## Arguments

- x:

  The return value of
  [`lingam_resit()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit.md)
  (a `ResitResult` object)

- threshold:

  Kept for interface consistency with the other
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) methods
  (default: 0); with 0/1 entries any value in `[0, 1)` returns all edges

- ...:

  Unused

## Value

data.frame(from, to, estimate)

## Examples

``` r
# \donttest{
if (requireNamespace("mgcv", quietly = TRUE)) {
  nonlinear <- generate_resit_sample(n = 300, seed = 1)
  model <- lingam_resit(nonlinear$data)
  tidy(model)
}
#>   from to estimate
#> 1   x0 x1        1
#> 2   x0 x2        1
#> 3   x1 x2        1
#> 4   x2 x3        1
# }
```
