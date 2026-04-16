
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DirectLiNGAM

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
<!-- badges: end -->

LiNGAM is a new method for estimating structural equation models or
linear Bayesian networks. It is based on using the non-Gaussianity of
the data.

This package is a port of the Python lingam package to R.

- [The LiNGAM Project](https://sites.google.com/view/sshimizu06/lingam)
- [lingam](https://github.com/cdt15/lingam)

`DirectLiNGAM` は、Pythonで公開されている [LiNGAM
パッケージ](https://github.com/cdt15/lingam) (LiNGAM: Linear
Non-Gaussian Acyclic Model) の R 言語への移植版です。

現在、開発中のアルファ版であり、動作確認とフィードバックを目的として公開しています。

## 特徴

- Direct LiNGAM アルゴリズムの実装
- ブートストラップ法による因果構造の安定性評価
- DiagrammeR を用いた推定結果の可視化

## Installation

You can install the development version of DirectLiNGAM from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("morimotoosamu/DirectLiNGAM")
```

## Requirements

- DiagrammeR

## Example

サンプルデータの呼び出しとDirect LiNGAMの実行。

``` r
library(DirectLiNGAM)
data(LiNGAM_sample)

model <- direct_lingam(LiNGAM_sample)
```

推定された因果順序。

``` r
cat("Index: ", model$causal_order, "\n")
#> Index:  4 3 1 6 5 8 2 9 7 10
cat("Names: ", colnames(LiNGAM_sample)[model$causal_order], "\n")
#> Names:  x3 x2 x0 x5 x4 x7 x1 x8 x6 x9
```

推定された隣接行列

``` r
B_hat <- model$adjacency_matrix
colnames(B_hat) <- rownames(B_hat) <- colnames(LiNGAM_sample)
cat("\n--- Estimated Adjacency Matrix ---\n")
#> 
#> --- Estimated Adjacency Matrix ---
print(round(B_hat, 3))
#>        x0     x1     x2     x3     x4     x5    x6     x7    x8 x9
#> x0  0.000  0.000 -0.010  3.053  0.000  0.000 0.000  0.000 0.000  0
#> x1  2.973  0.000  1.985  0.094 -0.005  0.002 0.000  0.004 0.000  0
#> x2  0.000  0.000  0.000  5.994  0.000  0.000 0.000  0.000 0.000  0
#> x3  0.000  0.000  0.000  0.000  0.000  0.000 0.000  0.000 0.000  0
#> x4  8.064  0.000 -1.001 -0.003  0.000 -0.016 0.000  0.000 0.000  0
#> x5  4.010  0.000 -0.018  0.080  0.000  0.000 0.000  0.000 0.000  0
#> x6 -0.068  1.983  0.024  0.050  0.013 -0.003 0.000 -0.003 0.010  0
#> x7  3.114  0.000 -0.025  0.083  1.487 -0.004 0.000  0.000 0.000  0
#> x8  0.047 -0.005  0.522 -0.088  0.011  2.011 0.000 -0.011 0.000  0
#> x9 -0.019  0.015 -0.017  7.067 -0.026 -0.036 0.993  0.016 0.015  0
```

推定された因果グラフの描画

``` r
plot_adjacency_diagrammer(
  B_hat,
  threshold = 1,
  labels = colnames(LiNGAM_sample),
  graph_label = "Estimated Causal Structure (Direct LiNGAM)",
  rankdir = "TB",
  shape = "ellipse",
  fillcolor = "lightyellow"
)
```

<img src="man/figures/README-plot_adjacency-1.png" alt="" width="100%" />

総合因果効果の推定

``` r
estimate_all_total_effects(LiNGAM_sample, model)
#>           x0          x1          x2        x3            x4           x5
#> x0  0.000000  0.00000000 -0.01006592  2.993077  0.0000000000  0.000000000
#> x1  3.000285  0.00000000  1.95406587 20.969321  0.0008783329  0.001802922
#> x2  0.000000  0.00000000  0.00000000  5.994427  0.0000000000  0.000000000
#> x3  0.000000  0.00000000  0.00000000  0.000000  0.0000000000  0.000000000
#> x4  7.999851  0.00000000 -1.08091865 17.943540  0.0000000000 -0.015895887
#> x5  4.010087  0.00000000 -0.05844897 11.973759  0.0000000000  0.000000000
#> x6  6.012803  1.98298503  3.89486773 41.944885  0.0104711199  0.020833822
#> x7 14.997303  0.00000000 -1.66407673 35.895240  1.4873512243 -0.027516599
#> x8  8.017238 -0.00543654  0.39927394 26.943909 -0.0052182213  2.010616291
#> x9  6.009163  1.98431849  3.89006150 48.971056  0.0080309484  0.015891490
#>           x6           x7         x8 x9
#> x0 0.0000000  0.000000000 0.00000000  0
#> x1 0.0000000  0.004071862 0.00000000  0
#> x2 0.0000000  0.000000000 0.00000000  0
#> x3 0.0000000  0.000000000 0.00000000  0
#> x4 0.0000000  0.000000000 0.00000000  0
#> x5 0.0000000  0.000000000 0.00000000  0
#> x6 0.0000000  0.005135646 0.01029517  0
#> x7 0.0000000  0.000000000 0.00000000  0
#> x8 0.0000000 -0.011047144 0.00000000  0
#> x9 0.9931672  0.021108335 0.02570678  0
```

事前知識を用いた推定

誤差独立性の p 値

## 注意事項

本パッケージは Python 版のすべての機能を網羅しているわけではありません。

本パッケージには Python 版に存在しない機能も追加されています。

## ライセンス

MIT License - オリジナルの作者である Shohei Shimizu
氏らに敬意を表します。

## フィードバック

バグの報告や機能のリクエストは、GitHub の Issues までお寄せください
