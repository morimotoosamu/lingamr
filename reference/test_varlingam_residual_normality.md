# Test the non-Gaussianity of VAR-LiNGAM residuals

LiNGAM assumes the error terms are non-Gaussian, so rejecting normality
(small p-value) supports the model assumption. By default the test is
run on the LiNGAM innovations `e_t = (I - B0) n_t` (the independent
errors the model assumes), where `n_t` are the stored VAR residuals; set
`on = "var"` to test the reduced-form VAR residuals `n_t` directly
instead.

## Usage

``` r
test_varlingam_residual_normality(
  result,
  method = "shapiro",
  alpha = 0.05,
  on = c("innovations", "var")
)
```

## Arguments

- result:

  a `VARLiNGAMResult` from
  [`lingam_var()`](https://morimotoosamu.github.io/lingamr/reference/lingam_var.md)

- method:

  normality test ("shapiro", "ks", "ad", "lillie", "jb"); see
  [`test_residual_normality()`](https://morimotoosamu.github.io/lingamr/reference/test_residual_normality.md)
  for package requirements

- alpha:

  significance level (default 0.05)

- on:

  which series to test: "innovations" (default, `e_t = (I - B0) n_t`) or
  "var" (the reduced-form VAR residuals `n_t`)

## Value

a `lingam_normality_test` data frame (one row per variable), printed via
[`print.lingam_normality_test()`](https://morimotoosamu.github.io/lingamr/reference/print.lingam_normality_test.md).

## References

Residual non-Gaussianity diagnostics inspired by the VARLiNGAM R code
(Gauss_Tests) of Moneta, A., Entner, D., Hoyer, P. O., & Coad, A.
(2013), *Oxford Bulletin of Economics and Statistics*, 75(5), 705-730.
<https://sites.google.com/site/dorisentner/publications/VARLiNGAM>

## Examples

``` r
s <- generate_varlingam_sample(n = 1000, seed = 42)
m <- lingam_var(s$data, lags = 1, reg_method = "ols", prune = FALSE)
test_varlingam_residual_normality(m)
#> === Residual Normality Test ===
#> Method:         shapiro
#> Sample size:    999
#> Significance:   0.050
#> Non-Gaussian:   3 / 3 variables
#> 
#>  variable statistic   p_value is_non_gauss skewness kurtosis
#>        x0    0.9498 < 2.2e-16         TRUE    0.088   -1.220
#>        x1    0.9536 < 2.2e-16         TRUE   -0.009   -1.238
#>        x2    0.9544 < 2.2e-16         TRUE   -0.045   -1.209
#> 
#> Interpretation:
#>   is_non_gauss = TRUE  -> rejects normality (supports LiNGAM assumption)
#>   is_non_gauss = FALSE -> cannot reject normality (LiNGAM may not fit)
#> 
#> All residuals are non-Gaussian. LiNGAM assumption is supported.
```
