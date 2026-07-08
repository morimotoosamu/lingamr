# Get a one-row summary of a LiMResult

Like
[`glance.LingamResult()`](https://morimotoosamu.github.io/lingamr/reference/glance.LingamResult.md),
with an additional `n_discrete` column giving the number of discrete
variables in the model.

## Usage

``` r
# S3 method for class 'LiMResult'
glance(x, ...)
```

## Arguments

- x:

  The return value of
  [`lingam_lim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_lim.md)
  (a `LiMResult` object)

- ...:

  Unused

## Value

A one-row data.frame(n_variables, n_edges, n_discrete, causal_order)

## Examples

``` r
set.seed(1)
dat <- generate_lim_sample(n = 300)
model <- lingam_lim(dat$data, is_continuous = dat$is_continuous)
glance(model)
#>   n_variables n_edges n_discrete   causal_order
#> 1           3       2          1 x1 -> x2 -> x3
```
