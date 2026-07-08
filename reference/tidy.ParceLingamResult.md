# Convert a ParceLingamResult to a tidy data.frame

Converts the estimated adjacency matrix into a long-format data.frame
with one edge per row, like
[`tidy.LingamResult()`](https://morimotoosamu.github.io/lingamr/reference/tidy.LingamResult.md).
`NA` entries of the adjacency matrix (variable pairs whose order could
not be resolved / suspected latent confounding) are kept as rows with
`estimate = NA` so they remain visible; drop them with e.g.
`subset(tidy(x), !is.na(estimate))` if not needed.

## Usage

``` r
# S3 method for class 'ParceLingamResult'
tidy(x, threshold = 0, ...)
```

## Arguments

- x:

  The return value of
  [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)
  (a `ParceLingamResult` object)

- threshold:

  Coefficients with an absolute value at or below this are not treated
  as edges (default: 0). `NA` entries are always kept.

- ...:

  Unused

## Value

data.frame(from, to, estimate)

## Examples

``` r
dat <- generate_parce_sample(n = 500, seed = 42)
model <- lingam_parce(dat$data)
tidy(model)
#>   from to   estimate
#> 1   x0 x1  0.5563819
#> 2   x0 x4  0.5680799
#> 3   x0 x5  0.5259416
#> 4   x2 x1  0.4714245
#> 5   x2 x3  0.7694629
#> 6   x2 x4 -0.5430309
#> 7   x3 x0  0.4913871
```
