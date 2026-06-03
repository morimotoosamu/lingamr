# 全変数間の総合因果効果を一括推定

全変数間の総合因果効果を一括推定

## Usage

``` r
estimate_all_total_effects(X, lingam_result, method = "lasso", lambda = "AIC")
```

## Arguments

- X:

  元データ (n_samples x n_features)

- lingam_result:

  lingam_direct() の返り値

- method:

  回帰手法 ("ols", "lasso", "adaptive_lasso")

- lambda:

  ラムダ選択 ("lambda.min", "lambda.1se", "AIC", "BIC")

## Value

総合因果効果の行列 (行: 結果変数, 列: 原因変数)

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

model <- LiNGAM_sample_1000$data |>
  lingam_direct()

LiNGAM_sample_1000$data |>
  estimate_all_total_effects(model)
#>          x0 x1        x2        x3       x4        x5
#> x0 0.000000  0  0.000000  3.033460 0.000000 0.0000000
#> x1 2.909144  0  1.889122 21.058733 0.000000 0.1550764
#> x2 0.000000  0  0.000000  5.992677 0.000000 0.0000000
#> x3 0.000000  0  0.000000  0.000000 0.000000 0.0000000
#> x4 8.001353  0 -1.308542 18.276121 0.000000 0.0000000
#> x5 4.014107  0  0.000000 12.179395 0.003005 0.0000000
```
