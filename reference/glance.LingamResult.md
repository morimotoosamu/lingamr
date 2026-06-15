# Get a one-row summary of a LingamResult

Summarizes the entire model in a single row. The data `X` is not
required because no residuals are computed. If residual-based
diagnostics are needed, use
[`summary_lingam()`](https://morimotoosamu.github.io/lingamr/reference/summary_lingam.md)
instead.

## Usage

``` r
# S3 method for class 'LingamResult'
glance(x, ...)
```

## Arguments

- x:

  The return value of
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
  (a `LingamResult` object)

- ...:

  Unused

## Value

A one-row data.frame(n_variables, n_edges, causal_order)

## Examples

``` r
dat <- generate_lingam_sample_6()
model <- lingam_direct(dat$data, reg_method = "ols")
glance(model)
#>   n_variables n_edges                     causal_order
#> 1           6      15 x3 -> x2 -> x0 -> x4 -> x5 -> x1
```
