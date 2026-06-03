# plot QQ

plot QQ

## Usage

``` r
plot_residual_qq(X, lingam_result, ncol = 3, nrow = NULL)
```

## Arguments

- X:

  元データ (matrix or data.frame)

- lingam_result:

  lingam_direct() の返り値

- ncol:

  Number of columns.

- nrow:

  Number of rows.

## Examples

``` r
# サンプルデータの呼び出し
LiNGAM_sample_1000 <- generate_lingam_sample_6()

# Direct LiNGAM の実行
result <- lingam_direct(LiNGAM_sample_1000$data)

plot_residual_qq(LiNGAM_sample_1000$data, result)
```
