# Get a one-row summary of a MultiGroupLingamResult

Summarizes the joint model in a single row. The causal order is shared
across groups; per-group edge counts are available via
`glance(get_group_result(x, i))`.

## Usage

``` r
# S3 method for class 'MultiGroupLingamResult'
glance(x, ...)
```

## Arguments

- x:

  The return value of
  [`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md)
  (a `MultiGroupLingamResult` object)

- ...:

  Unused

## Value

A one-row data.frame(n_groups, n_variables, causal_order)

## Examples

``` r
mg <- generate_multi_group_sample()
model <- lingam_multi_group(mg$data_list, reg_method = "ols")
glance(model)
#>   n_groups n_variables                     causal_order
#> 1        2           6 x3 -> x0 -> x5 -> x2 -> x4 -> x1
```
