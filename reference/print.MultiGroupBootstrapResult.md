# Print method for MultiGroupBootstrapResult

Print method for MultiGroupBootstrapResult

## Usage

``` r
# S3 method for class 'MultiGroupBootstrapResult'
print(x, ...)
```

## Arguments

- x:

  MultiGroupBootstrapResult object

- ...:

  Additional arguments (unused)

## Value

The input object `x`, invisibly.

## Examples

``` r
mg <- generate_multi_group_sample()
bs <- lingam_multi_group_bootstrap(mg$data_list,
  n_sampling = 10L, reg_method = "ols", seed = 42
)
#> Multi-group bootstrap: 10 iterations, 2 groups, method=ols (sequential)
#>   iteration 1 / 10
#>   iteration 10 / 10
#> Completed in 0.2 seconds.
print(bs)
#> MultiGroupBootstrapResult: 2 groups
#>   [group1] 10 samplings, 6 features
#>   [group2] 10 samplings, 6 features
```
