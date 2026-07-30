# Run several normality tests on VARMA-LiNGAM residuals at once

Convenience wrapper (analogous to the Moneta `Gauss_Tests`) that applies
multiple normality tests to the residuals and returns a single table
with one p-value column per method plus per-variable skewness and excess
kurtosis. Methods whose optional package is unavailable are skipped with
a warning.

## Usage

``` r
test_varmalingam_residual_normality_all(
  result,
  methods = c("shapiro", "ad", "lillie", "jb"),
  alpha = 0.05,
  on = c("innovations", "varma")
)
```

## Arguments

- result:

  a `VARMALiNGAMResult` from
  [`lingam_varma()`](https://morimotoosamu.github.io/lingamr/reference/lingam_varma.md)

- methods:

  character vector of tests to run; any of "shapiro", "ks", "ad",
  "lillie", "jb" (default runs shapiro/ad/lillie/jb)

- alpha:

  significance level (default 0.05)

- on:

  which series to test: "innovations" (default) or "varma"

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
s <- generate_varmalingam_sample(n = 1000, seed = 42)
m <- lingam_varma(s$data,
  order = c(1, 1), criterion = NULL,
  reg_method = "ols", prune = FALSE
)
test_varmalingam_residual_normality_all(m, methods = c("shapiro", "jb"))
#>   variable     skewness  kurtosis    p_shapiro         p_jb all_non_gauss
#> 1       x0  0.087734895 -1.221621 5.049812e-18 1.709743e-14          TRUE
#> 2       x1  0.005493456 -1.229598 6.460986e-17 2.153833e-14          TRUE
#> 3       x2 -0.042210104 -1.206647 5.256731e-17 5.961898e-14          TRUE
```
