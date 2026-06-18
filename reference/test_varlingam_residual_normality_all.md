# Run several normality tests on VAR-LiNGAM residuals at once

Convenience wrapper (analogous to the Moneta `Gauss_Tests`) that applies
multiple normality tests to the residuals and returns a single table
with one p-value column per method plus per-variable skewness and excess
kurtosis. Methods whose optional package is unavailable are skipped with
a warning.

## Usage

``` r
test_varlingam_residual_normality_all(
  result,
  methods = c("shapiro", "ad", "lillie", "jb"),
  alpha = 0.05,
  on = c("innovations", "var")
)
```

## Arguments

- result:

  a `VARLiNGAMResult` from
  [`lingam_var()`](https://morimotoosamu.github.io/lingamr/reference/lingam_var.md)

- methods:

  character vector of tests to run; any of "shapiro", "ks", "ad",
  "lillie", "jb" (default runs shapiro/ad/lillie/jb)

- alpha:

  significance level (default 0.05)

- on:

  which series to test: "innovations" (default) or "var"

## Value

a data frame with columns `variable`, `skewness`, `kurtosis`, one
`p_<method>` column per method, and `all_non_gauss` (TRUE when every run
test rejects normality for that variable).

## References

Analogous to the multi-test residual check (Gauss_Tests) in the
VARLiNGAM R code of Moneta, A., Entner, D., Hoyer, P. O., & Coad, A.
(2013), *Oxford Bulletin of Economics and Statistics*, 75(5), 705-730.
<https://sites.google.com/site/dorisentner/publications/VARLiNGAM>

## Examples

``` r
s <- generate_varlingam_sample(n = 1000, seed = 42)
m <- lingam_var(s$data, lags = 1, reg_method = "ols", prune = FALSE)
test_varlingam_residual_normality_all(m, methods = c("shapiro", "jb"))
#> Registered S3 method overwritten by 'quantmod':
#>   method            from
#>   as.zoo.data.frame zoo 
#>   variable     skewness  kurtosis    p_shapiro         p_jb all_non_gauss
#> 1       x0  0.088013433 -1.219504 6.041150e-18 1.898481e-14          TRUE
#> 2       x1 -0.008502543 -1.237896 3.173909e-17 1.398881e-14          TRUE
#> 3       x2 -0.045078629 -1.209086 4.530246e-17 5.162537e-14          TRUE
```
