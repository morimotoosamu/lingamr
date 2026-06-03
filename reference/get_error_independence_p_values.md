# 誤差の独立性検定の p 値を計算

誤差の独立性検定の p 値を計算

## Usage

``` r
get_error_independence_p_values(X, lingam_result, method = "spearman")
```

## Arguments

- X:

  元データ (matrix or data.frame)

- lingam_result:

  lingam_direct() の返り値

- method:

  相関係数の種類 ("spearman", "pearson", "kendall")

## Value

p 値の行列 (n_features x n_features)

## Examples

``` r
# サンプルデータの呼び出し
LiNGAM_sample_1000 <- generate_lingam_sample_6()

# Direct LiNGAM の実行
result <- LiNGAM_sample_1000$data |>
  lingam_direct()

# p 値の計算（デフォルト: Spearman）
p_vals <- get_error_independence_p_values(LiNGAM_sample_1000$data, result)
round(p_vals, 3)
#>       x0    x1    x2    x3    x4    x5
#> x0    NA 0.988 0.214 0.976 0.484 0.954
#> x1 0.988    NA 0.986 0.991 0.323 0.882
#> x2 0.214 0.986    NA 0.919 0.100 0.124
#> x3 0.976 0.991 0.919    NA 0.806 0.974
#> x4 0.484 0.323 0.100 0.806    NA 0.643
#> x5 0.954 0.882 0.124 0.974 0.643    NA

# Kendall で計算
p_vals_k <- get_error_independence_p_values(LiNGAM_sample_1000$data, result, method = "kendall")
round(p_vals_k, 3)
#>       x0    x1    x2    x3    x4    x5
#> x0    NA 0.986 0.225 0.996 0.478 0.961
#> x1 0.986    NA 0.978 0.969 0.320 0.894
#> x2 0.225 0.978    NA 0.912 0.104 0.131
#> x3 0.996 0.969 0.912    NA 0.798 0.954
#> x4 0.478 0.320 0.104 0.798    NA 0.641
#> x5 0.961 0.894 0.131 0.954 0.641    NA
```
