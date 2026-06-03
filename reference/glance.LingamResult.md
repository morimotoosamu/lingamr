# LingamResult の1行サマリを取得

モデル全体を1行に要約する。残差を計算しないためデータ `X` は不要。
残差ベースの診断が必要な場合は
[`summary_lingam()`](https://morimotoosamu.github.io/lingamr/reference/summary_lingam.md)
を使用すること。

## Usage

``` r
# S3 method for class 'LingamResult'
glance(x, ...)
```

## Arguments

- x:

  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
  の返り値（`LingamResult` オブジェクト）

- ...:

  未使用

## Value

1 行の data.frame(n_variables, n_edges, causal_order)

## Examples

``` r
dat <- generate_lingam_sample_6()
model <- lingam_direct(dat$data, reg_method = "ols")
glance(model)
#>   n_variables n_edges                     causal_order
#> 1           6      15 x3 -> x2 -> x0 -> x4 -> x5 -> x1
```
