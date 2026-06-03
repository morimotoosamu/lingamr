# 全変数間の総合因果効果を一括推定

全変数間の総合因果効果を一括推定

## Usage

``` r
estimate_all_total_effects(
  X,
  lingam_result,
  method = "adaptive_lasso",
  lambda = "BIC"
)
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

総合因果効果の行列 (n_features x n_features)。 **規則: `TE[i, j]` は変数
j から変数 i への総合因果効果（j → i）。** 隣接行列 `adjacency_matrix`
と同じ添字規則。直接効果と間接効果の合計。

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

model <- LiNGAM_sample_1000$data |>
  lingam_direct()

LiNGAM_sample_1000$data |>
  estimate_all_total_effects(model)
#>          x0 x1        x2        x3 x4 x5
#> x0 0.000000  0  0.000000  3.033460  0  0
#> x1 2.896907  0  1.909712 21.058733  0  0
#> x2 0.000000  0  0.000000  5.992677  0  0
#> x3 0.000000  0  0.000000  0.000000  0  0
#> x4 8.001464  0 -1.308131 18.276121  0  0
#> x5 4.015103  0  0.000000 12.179395  0  0
```
