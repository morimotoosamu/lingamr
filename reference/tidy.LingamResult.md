# LingamResult を tidy な data.frame に変換

推定された隣接行列を、1 行が 1 エッジの long 形式 data.frame
に変換する。 `B[i, j]`（j → i の係数）の規則に従い、`from`
列が原因、`to` 列が結果となる。 ggplot2 や ggraph での可視化、dplyr
でのフィルタリングに便利。

## Usage

``` r
# S3 method for class 'LingamResult'
tidy(x, threshold = 0, ...)
```

## Arguments

- x:

  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
  の返り値（`LingamResult` オブジェクト）

- threshold:

  この絶対値以下の係数はエッジとみなさない (default: 0)

- ...:

  未使用

## Value

data.frame(from, to, estimate)。`from`/`to` は変数名（文字列）、
`estimate` は因果係数。エッジが無ければ 0 行の data.frame。

## Examples

``` r
dat <- generate_lingam_sample_6()
model <- lingam_direct(dat$data, reg_method = "ols")
tidy(model)
#>    from to     estimate
#> 1    x0 x1  3.236756411
#> 2    x0 x4  7.992316238
#> 3    x0 x5  3.873373532
#> 4    x2 x0 -0.040123987
#> 5    x2 x1  1.965485693
#> 6    x2 x4 -1.061625935
#> 7    x2 x5  0.069075201
#> 8    x3 x0  3.273910192
#> 9    x3 x1  0.013952441
#> 10   x3 x2  5.992677091
#> 11   x3 x4  0.393730192
#> 12   x3 x5 -0.314606489
#> 13   x4 x1 -0.033970558
#> 14   x4 x5  0.018074531
#> 15   x5 x1  0.005627409
```
