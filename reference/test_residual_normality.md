# Test normality of residuals from Direct LiNGAM

Calculate residuals (error terms) from the estimated adjacency matrix
and test their normality. Since LiNGAM assumes non-Gaussian errors,
rejecting normality (small p-value) supports the LiNGAM model
assumption.

## Usage

``` r
test_residual_normality(X, lingam_result, method = "shapiro", alpha = 0.05)
```

## Arguments

- X:

  original data matrix or data.frame

- lingam_result:

  result from lingam_direct()

- method:

  normality test method "shapiro" : Shapiro-Wilk test (default, n
  \<= 5000) "ks" : Kolmogorov-Smirnov test (n \> 5000) "ad" :
  Anderson-Darling test (requires nortest package) "lillie" : Lilliefors
  test (requires nortest package) "jb" : Jarque-Bera test (requires
  tseries package)

- alpha:

  significance level (default: 0.05)

## Value

data.frame with test results for each variable

## Examples

``` r
# サンプルデータの呼び出し
LiNGAM_sample_1000 <- generate_lingam_sample_6()

# Direct LiNGAM の実行
result <- lingam_direct(LiNGAM_sample_1000$data)

# Shapiro-Wilk (default)
test_residual_normality(LiNGAM_sample_1000$data, result)
#> === Residual Normality Test ===
#> Method:         shapiro
#> Sample size:    1000
#> Significance:   0.050
#> Non-Gaussian:   6 / 6 variables
#> 
#>  variable statistic   p_value is_non_gauss skewness kurtosis
#>        x0    0.9516 < 2.2e-16         TRUE    0.061   -1.215
#>        x1    0.9521 < 2.2e-16         TRUE    0.026   -1.213
#>        x2    0.9557 < 2.2e-16         TRUE    0.083   -1.170
#>        x3    0.9578  2.25e-16         TRUE    0.025   -1.163
#>        x4    0.9546 < 2.2e-16         TRUE   -0.003   -1.205
#>        x5    0.9536 < 2.2e-16         TRUE   -0.052   -1.206
#> 
#> Interpretation:
#>   is_non_gauss = TRUE  -> rejects normality (supports LiNGAM assumption)
#>   is_non_gauss = FALSE -> cannot reject normality (LiNGAM may not fit)
#> 
#> All residuals are non-Gaussian. LiNGAM assumption is supported.
```
