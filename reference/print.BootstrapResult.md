# Display the contents of a BootstrapResult

Display the contents of a BootstrapResult

## Usage

``` r
# S3 method for class 'BootstrapResult'
print(x, ...)
```

## Arguments

- x:

  BootstrapResult object

- ...:

  Additional arguments (for S3 method compatibility)

## Value

The input object `x`, invisibly.

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

bs_model <- lingam_direct_bootstrap(LiNGAM_sample_1000$data, n_sampling = 30L, seed = 42)
#> Bootstrap: 30 iterations, method=adaptive_lasso (sequential)
#>   iteration 1 / 30
#>   iteration 10 / 30
#>   iteration 20 / 30
#>   iteration 30 / 30
#> Completed in 0.9 seconds.

print(bs_model)
#> BootstrapResult: 30 samplings, 6 features
```
