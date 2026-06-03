# 指定した2変数間の総合因果効果を推定

指定した2変数間の総合因果効果を推定

## Usage

``` r
estimate_total_effect(
  X,
  lingam_result,
  from_index,
  to_index,
  method = "adaptive_lasso",
  lambda = "BIC"
)
```

## Arguments

- X:

  元データ (matrix or data.frame)

- lingam_result:

  lingam_direct() の返り値

- from_index:

  原因変数 (1-based index or 変数名)

- to_index:

  結果変数 (1-based index or 変数名)

- method:

  回帰手法 ("ols", "lasso", "adaptive_lasso")デフォルトはadaptive_lasso

- lambda:

  ラムダ選択 ("lambda.min", "lambda.1se", "AIC", "BIC",
  "oracle")デフォルトはBIC

## Value

推定された総合因果効果

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

model <- LiNGAM_sample_1000$data |>
  lingam_direct()

LiNGAM_sample_1000$data |>
  estimate_total_effect(model, 4, 1)
#>      x3 
#> 3.03346 
```
