# lingamr

[![R-CMD-check](https://github.com/morimotoosamu/lingamr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/morimotoosamu/lingamr/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/lingamr)](https://CRAN.R-project.org/package=lingamr)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)

LiNGAMは、データの非ガウス性を利用して構造方程式モデルや線形ベイジアンネットワークを推定する手法です。

`lingamr` は、Pythonで公開されている
[LiNGAM](https://github.com/cdt15/lingam) パッケージ(LiNGAM: Linear
Non-Gaussian Acyclic Model)をRに移植したものです。

English version:
[README.md](https://github.com/morimotoosamu/lingamr/blob/HEAD/README.md)

- [The LiNGAM Project](https://sites.google.com/view/sshimizu06/lingam)
- [lingam (Python)](https://github.com/cdt15/lingam)

## 特徴

- Direct
  LiNGAM([`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md))。隣接行列の推定に複数の回帰バックエンド
  (OLS、LASSO、適応的LASSO、ridge)を選択可能
- VAR-LiNGAM([`lingam_var()`](https://morimotoosamu.github.io/lingamr/reference/lingam_var.md))による時系列データの因果探索
- 潜在交絡変数に頑健なアルゴリズム: BottomUpParceLiNGAM
  ([`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md))とRCD([`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md))
- 複数データセットで共通の因果順序を同時推定するMultiGroupDirectLiNGAM
  ([`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md))
- 高次元データ(`p`が大きい、または`p > n`)向けのHighDimDirectLiNGAM
  ([`lingam_high_dim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_high_dim.md))
- 連続・二値の混合データに対応するLiM([`lingam_lim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_lim.md))
- 欠測値を含むデータの因果探索を行う[`bootstrap_with_imputation()`](https://morimotoosamu.github.io/lingamr/reference/bootstrap_with_imputation.md)と、
  lavaanによるSEM適合度評価[`evaluate_model_fit()`](https://morimotoosamu.github.io/lingamr/reference/evaluate_model_fit.md)
- ブートストラップ法による因果構造の安定性評価(因果順序の安定性を含む)
- モデル診断:
  残差の独立性・正規性検定と、一括実行できる[`summary_lingam()`](https://morimotoosamu.github.io/lingamr/reference/summary_lingam.md)
- DiagrammeR(インタラクティブ)とggplot2の`autoplot()`(静的)による可視化
- broomスタイルのtidier([`tidy()`](https://generics.r-lib.org/reference/tidy.html)
  / [`glance()`](https://generics.r-lib.org/reference/glance.html))

Python版のすべての機能を含んでいるわけではなく、逆にPython版にはない機能も
含まれています。

## インストール

CRANから`lingamr`をインストールできます。

``` r

install.packages("lingamr")
```

または、[GitHub](https://github.com/morimotoosamu/lingamr)から開発版を
インストールできます。

``` r

# install.packages("pak")
pak::pak("morimotoosamu/lingamr")
```

一部の機能は次のSuggestsパッケージに依存します: `DiagrammeR`
(インタラクティブなプロット)、`igraph`と`ggplot2`(静的な`autoplot()`
グラフとQQプロット)、`glmnet`(適応的LASSO)、`nortest` / `tseries`
(残差検定)。

## クイックスタート

``` r

library(lingamr)

# 6変数のLiNGAMモデルからサンプルデータを生成
x <- generate_lingam_sample_6(n = 1000)

# Direct LiNGAMで因果構造を推定
model <- lingam_direct(x$data)

# 推定された因果順序(変数名で表示)
colnames(x$data)[model$causal_order]
#> [1] "x3" "x2" "x0" "x4" "x5" "x1"
```

``` r

# 推定された因果グラフを可視化
model$adjacency_matrix |>
  plot_adjacency(
    labels    = colnames(model$adjacency_matrix),
    title     = "Estimated Causal Structure (Direct LiNGAM)",
    rankdir   = "TB",
    shape     = "ellipse",
    fillcolor = "lightgreen"
  )
```

![](reference/figures/README-quickstart_plot-1.png)

## 詳しく知るには

事前知識の組み込み、総因果効果、残差の独立性・正規性検定、ブートストラップ
(並列実行を含む)まで一通り解説したチュートリアルは、vignetteを参照して
ください。

``` r

vignette("lingamr")    # 英語版
vignette("lingamr-ja") # 日本語版
```

## ライセンス

MIT License

Original work: Copyright (c) 2019 T.Ikeuchi, G.Haraoka, M.Ide,
W.Kurebayashi, S.Shimizu

Portions of this work: Copyright (c) 2026 O.Morimoto

## 参考文献

### アルゴリズム

- Shimizu, S. et al. (2011). DirectLiNGAM: A direct method for learning
  a linear non-Gaussian structural equation model. *Journal of Machine
  Learning Research*, 12, 1225-1248.

### 元となったPython実装

- Ikeuchi, T. et al. (2023). Python package for LiNGAM algorithms.
  *Journal of Machine Learning Research*, 24(14), 1-7.
  <https://github.com/cdt15/lingam>

### 書籍

- 清水昌平(2017)『統計的因果探索』講談社.
- 梅津佑太・西村龍映・上田勇祐(2020)『スパース回帰分析とパターン認識』講談社.
- 鈴木譲(2025)『グラフィカルモデルと因果探索100問 with R』共立出版.

### 参考にしたRパッケージ

- G. Kikuchi (2020). rlingam <https://github.com/gkikuchi/rlingam>

## 謝辞

本パッケージの開発にあたり、以下の方々および団体にご支援いただきました。

- Hirohiko Asano
- Yoshiyuki Okuse
- Takahiko Umeyama
- JMRA(日本マーケティング・リサーチ協会)

本パッケージの開発にあたっては、AIコーディングツール(Google Geminiおよび
Anthropic Claude)の支援を受けました。生成されたコードはすべて著者が
レビュー・テスト・検証しています。

## フィードバック

バグ報告・機能要望はGitHub Issuesからお願いします。
