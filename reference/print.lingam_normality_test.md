# Print method for lingam_normality_test

Print method for lingam_normality_test

## Usage

``` r
# S3 method for class 'lingam_normality_test'
print(x, ...)
```

## Arguments

- x:

  lingam_normality_test object

- ...:

  additional arguments

## Value

The input object `x`, invisibly.

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()
result <- lingam_direct(LiNGAM_sample_1000$data, reg_method = "ols")
print(test_residual_normality(LiNGAM_sample_1000$data, result))
#> === Residual Normality Test ===
#> Method:         shapiro
#> Sample size:    1000
#> Significance:   0.050
#> Non-Gaussian:   6 / 6 variables
#> 
#>  variable statistic   p_value is_non_gauss skewness kurtosis
#>        x0    0.9526 < 2.2e-16         TRUE    0.063   -1.214
#>        x1    0.9533 < 2.2e-16         TRUE    0.026   -1.208
#>        x2    0.9557 < 2.2e-16         TRUE    0.083   -1.170
#>        x3    0.9578  2.25e-16         TRUE    0.025   -1.163
#>        x4    0.9556 < 2.2e-16         TRUE   -0.002   -1.205
#>        x5    0.9552 < 2.2e-16         TRUE   -0.044   -1.196
#> 
#> Interpretation:
#>   is_non_gauss = TRUE  -> rejects normality (supports LiNGAM assumption)
#>   is_non_gauss = FALSE -> cannot reject normality (LiNGAM may not fit)
#> 
#> All residuals are non-Gaussian. LiNGAM assumption is supported.
```
