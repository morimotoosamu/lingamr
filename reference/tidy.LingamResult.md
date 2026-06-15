# Convert a LingamResult to a tidy data.frame

Converts the estimated adjacency matrix into a long-format data.frame
with one edge per row. Following the `B[i, j]` convention (the
coefficient for j -\> i), the `from` column is the cause and the `to`
column is the effect. Convenient for visualization with ggplot2 or
ggraph and for filtering with dplyr.

## Usage

``` r
# S3 method for class 'LingamResult'
tidy(x, threshold = 0, ...)
```

## Arguments

- x:

  The return value of
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
  (a `LingamResult` object)

- threshold:

  Coefficients with an absolute value at or below this are not treated
  as edges (default: 0)

- ...:

  Unused

## Value

data.frame(from, to, estimate). `from`/`to` are variable names (strings)
and `estimate` is the causal coefficient. Returns a 0-row data.frame if
there are no edges.

## Examples

``` r
dat <- generate_lingam_sample_6()
model <- lingam_direct(dat$data, reg_method = "ols")
tidy(model)
#>    from to     estimate
#> 1    x0 x1  3.236756411
#> 2    x0 x4  7.992316238
#> 3    x0 x5  3.873373532
#> 4    x2 x0 -0.040123987
#> 5    x2 x1  1.965485693
#> 6    x2 x4 -1.061625935
#> 7    x2 x5  0.069075201
#> 8    x3 x0  3.273910192
#> 9    x3 x1  0.013952441
#> 10   x3 x2  5.992677091
#> 11   x3 x4  0.393730192
#> 12   x3 x5 -0.314606489
#> 13   x4 x1 -0.033970558
#> 14   x4 x5  0.018074531
#> 15   x5 x1  0.005627409
```
