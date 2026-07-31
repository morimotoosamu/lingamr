# 特殊なデータ: 混合変数・複数グループ・欠測値

Direct
LiNGAMは、連続変数からなる単一の完全なデータ行列を仮定する。この記事では、
その仮定が崩れるデータのための3つの拡張を扱う。

- **LiM**（[`lingam_lim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_lim.md)）:
  連続変数と離散変数（二値またはポアソンカウント）が 混在するデータ。
- **MultiGroup Direct
  LiNGAM**（[`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md)）:
  因果構造を共有するが 係数値は共有しない複数のデータセット。
- **欠測データ**（[`bootstrap_with_imputation()`](https://morimotoosamu.github.io/lingamr/reference/bootstrap_with_imputation.md)）:
  `NA` を含むデータのための ブートストラップ + 多重代入。

``` r

library(lingamr)
```

## 混合データのためのLiNGAM（LiM）

Direct
LiNGAMはすべての変数が連続であることを仮定する。[`lingam_lim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_lim.md)
はこの仮定を 緩和し、Zeng et al. (2022)
に従って、連続変数と離散変数が混在するデータ
から因果構造を推定する。NOTEARS流の連続最適化（「global」フェーズ）と、エッジの
方向・枝刈り・エッジ追加に関する組合せ的な局所探索（「local」フェーズ）を組み合わせ
たものである。

[`generate_lim_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_lim_sample.md)
は、連続変数と離散変数からなる既知の因果連鎖`x1`（連続） -\>
`x2`（離散）-\> `x3`（連続）を持つ小規模なデータセットを生成する。

``` r

set.seed(1)
lim_dat <- generate_lim_sample(n = 2000)
head(lim_dat$data)
#>           x1 x2         x3
#> 1  0.1182559  0 -0.8936636
#> 2 -1.8695490  0 -1.2651618
#> 3 -3.3867259  0  2.1530815
#> 4 -0.4395899  1  0.6618645
#> 5  0.3215812  0 -0.4652433
#> 6  1.6721779  0  0.8429619
lim_dat$is_continuous
#> [1]  TRUE FALSE  TRUE
```

[`lingam_lim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_lim.md)
には、各列が連続（`TRUE`）か離散（`FALSE`）かを示す論理ベクトル
`is_continuous`
が必要である。最適化はランダムな初期点から始まるため、再現性を持たせる
には [`set.seed()`](https://rdrr.io/r/base/Random.html) が必要である。

``` r

lim_result <- lingam_lim(lim_dat$data, is_continuous = lim_dat$is_continuous)
print(lim_result)
#> LiM Result
#>   Variables : 3
#>   Variable types: continuous, discrete, continuous
#>   Causal order: x1 -> x2 -> x3
#> 
#> Adjacency matrix (row = to, col = from):
#>    x1    x2 x3
#> x1  0 0.000  0
#> x2  1 0.000  0
#> x3  0 1.657  0
```

[`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
と同様、`adjacency_matrix` は `B[i, j]` = j -\> i という規約 （行 =
to、列 = from）に従い、`causal_order`
は推定されたトポロジカル順序を1始まり のインデックスで並べたものである。

``` r

colnames(lim_dat$data)[lim_result$causal_order]
#> [1] "x1" "x2" "x3"
```

### カウント変数（ポアソン）

離散変数は既定では二値（0/1）である。`is_poisson = TRUE`
を指定するとポアソン分布に
従うカウント変数として扱われる（localフェーズがポアソン回帰の対数尤度でスコアリング
する）。`generate_lim_sample(is_poisson = TRUE)`
で対応するカウントデータの例を生成 できる。

``` r

set.seed(2)
lim_pois <- generate_lim_sample(n = 2000, is_poisson = TRUE)
head(lim_pois$data)
#>            x1 x2         x3
#> 1  0.95924714  2  1.9627734
#> 2 -0.27872569  0 -0.3284517
#> 3 -1.63830626  0  0.4878718
#> 4 -0.20204595  1  1.7284788
#> 5  0.02041924  5  3.2550128
#> 6 -0.42685945  1  0.5478110

lim_pois_result <- lingam_lim(lim_pois$data,
  is_continuous = lim_pois$is_continuous, is_poisson = TRUE
)
print(lim_pois_result)
#> LiM Result
#>   Variables : 3
#>   Variable types: continuous, discrete (count), continuous
#>   Causal order: x1 -> x2 -> x3
#> 
#> Adjacency matrix (row = to, col = from):
#>    x1    x2 x3
#> x1  0 0.000  0
#> x2  1 0.000  0
#> x3  0 0.505  0
```

localフェーズのエッジ重みの規約や、Python実装との数値的な違いの詳細は
[`?lingam_lim`](https://morimotoosamu.github.io/lingamr/reference/lingam_lim.md)
を参照のこと。

## MultiGroup Direct LiNGAM

[`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
は単一のデータセットに適用するものである。データが複数のソース
から得られ、同じ因果構造を共有しているが効果の強さは同じではないと考えられる場合
（例えば、同じ研究を複数拠点で実施した場合や、同じプロセスを異なる時期に観測した
場合）、[`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md)
はShimizu (2012) に従い、各グループに固有の隣接行列
（構造係数）を許容しつつ、すべてのグループにわたって**共通の因果順序**を同時に推定
する。

[`generate_multi_group_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_multi_group_sample.md)
は、[`generate_lingam_sample_6()`](https://morimotoosamu.github.io/lingamr/reference/generate_lingam_sample_6.md)
の因果構造を共有
しつつ、グループごとに係数がわずかに異なる2つのデータセットを生成する。

``` r

mg <- generate_multi_group_sample(n = c(1000, 1000), seed = 42)
lapply(mg$data_list, head, 3)
#> $group1
#>         x0        x1       x2        x3        x4        x5
#> 1 2.814924 18.017120 4.543655 0.6333728 18.160090 12.236660
#> 2 1.889685 10.956005 2.188091 0.3175366 13.172754  7.932657
#> 3 1.008905  6.990652 1.953131 0.2409218  6.702107  4.797122
#> 
#> $group2
#>          x0        x1        x2         x3        x4        x5
#> 1 0.7259014  5.482225 0.9301592 0.01259095  5.275903  3.321061
#> 2 2.4321051 17.252303 3.3989705 0.41696287 16.459975 11.368701
#> 3 1.5550457 10.342355 1.8713591 0.24518297 10.520165  7.859683
```

``` r

mg_result <- lingam_multi_group(mg$data_list, reg_method = "ols")
print(mg_result)
#> Multi-Group Direct LiNGAM Result
#>   Groups      : 2 (group1, group2)
#>   Variables   : 6
#>   Causal order (common): x3 -> x0 -> x5 -> x2 -> x4 -> x1
#> 
#> [group1] Adjacency matrix (row = to, col = from):
#>        x0 x1     x2     x3     x4    x5
#> x0  0.000  0  0.000  3.033  0.000 0.000
#> x1  3.237  0  1.965  0.014 -0.034 0.006
#> x2 -0.236  0  0.000  6.112  0.000 0.049
#> x3  0.000  0  0.000  0.000  0.000 0.000
#> x4  7.921  0 -1.063  0.399  0.000 0.018
#> x5  4.016  0  0.000 -0.003  0.000 0.000
#> 
#> [group2] Adjacency matrix (row = to, col = from):
#>       x0 x1     x2     x3    x4     x5
#> x0 0.000  0  0.000  3.504 0.000  0.000
#> x1 2.732  0  2.568  0.083 0.034  0.093
#> x2 0.154  0  0.000  6.322 0.000 -0.024
#> x3 0.000  0  0.000  0.000 0.000  0.000
#> x4 8.483  0 -1.487 -0.110 0.000  0.006
#> x5 4.515  0  0.000 -0.045 0.000  0.000
```

`causal_order` はすべてのグループで共有される。`adjacency_matrices`
はグループ ごとに1つの行列を保持し、それぞれ通常の `B[i, j]` = j -\> i
という規約に従う。

`lingamr`
のシングルグループ向けツール群（総因果効果、独立性検定、プロット）で
単一グループを分析するには、[`get_group_result()`](https://morimotoosamu.github.io/lingamr/reference/get_group_result.md)
で通常の `LingamResult` として 取り出す。

``` r

g1 <- get_group_result(mg_result, "group1")
class(g1)
#> [1] "LingamResult"

estimate_all_total_effects(mg$data_list$group1, g1, method = "ols")
#>             x0 x1        x2        x3          x4          x5
#> x0  0.00000000  0  0.000000  3.033460  0.00000000  0.00000000
#> x1  2.90911952  0  2.001580 21.058733 -0.03397056  0.10299386
#> x2 -0.03933572  0  0.000000  5.992677  0.00000000  0.04894766
#> x3  0.00000000  0  0.000000  0.000000  0.00000000  0.00000000
#> x4  8.03407606  0 -1.062516 18.276121  0.00000000 -0.03416285
#> x5  4.01586857  0  0.000000 12.179395  0.00000000  0.00000000
```

[`lingam_multi_group_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group_bootstrap.md)
も、同じ結合的な方法でブートストラップの安定性
推定を提供する。各反復で各グループを個別にリサンプルした上で、因果順序とグループ
ごとの隣接行列を結合的に再推定する。グループごとの `BootstrapResult`
オブジェクト
からなる名前付きリストを返すため、既存のブートストラップ照会関数をグループごとに
そのまま適用できる。

``` r

mg_bs <- lingam_multi_group_bootstrap(mg$data_list,
  n_sampling = 20L, reg_method = "ols", seed = 1, verbose = FALSE
)
get_probabilities(mg_bs$group1)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]  0.0  0.0 0.30    1 0.00 0.00
#> [2,]  1.0  0.0 1.00    1 0.70 0.80
#> [3,]  0.7  0.0 0.00    1 0.00 0.45
#> [4,]  0.0  0.0 0.00    0 0.00 0.00
#> [5,]  1.0  0.3 1.00    1 0.00 0.55
#> [6,]  1.0  0.2 0.55    1 0.45 0.00
```

[`lingam_multi_group_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group_bootstrap.md)
の総因果効果は、回帰によってではなく、各反復の
隣接行列に対するパス係数の積として計算される点に注意する。これは上流のPython実装
と一致するが、[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md)
の回帰ベースの
[`estimate_all_total_effects()`](https://morimotoosamu.github.io/lingamr/reference/estimate_all_total_effects.md)
とは異なる。

## 欠測データによる因果探索

これまでのアルゴリズムはすべて、完全なデータ行列を仮定している。`X`
に欠測値
（`NA`）が含まれる場合、[`bootstrap_with_imputation()`](https://morimotoosamu.github.io/lingamr/reference/bootstrap_with_imputation.md)
はブートストラップリサンプ
リングと多重代入を組み合わせる。各リサンプルは複数の完全なデータセットに補完され、
[`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md)（補完されたコピーを、1つの因果順序を共有する「グループ」と
して扱う）を用いて、それらにわたって共通の因果構造が結合的に推定される。これは
Pythonの `lingam.tools.bootstrap_with_imputation()` のR移植版である。

``` r

sample6_na <- generate_lingam_sample_6(n = 1000, seed = 1)
X_na <- sample6_na$data
set.seed(1)
X_na$x5[sample.int(nrow(X_na), size = round(0.1 * nrow(X_na)))] <- NA # MCAR 10% on x5
```

``` r

bwi <- bootstrap_with_imputation(X_na,
  n_sampling = 20L, n_repeats = 5L, seed = 42, verbose = FALSE
)
print(bwi)
#> ImputationBootstrapResult: 20 samplings x 5 repeats, 6 features, 100 missing cells (original data)
```

デフォルトの補完器は
`mice::mice(method = "norm")`（ベイズ線形回帰）であり、上流の
Pythonのデフォルト（`IterativeImputer(sample_posterior = TRUE)`）に最も近い標準的な
R実装である。数値結果はPython実装とは一致しない。補完器と因果探索ステップの両方
とも、`imputer` 引数と `cd_fit` 引数を通じてカスタム `function`
に差し替えることが できる。

各反復で（補完データセットごとに1つずつ）`n_repeats`
個の隣接行列が生成されるため、 結果の形状は
[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md)
とは異なる。[`as_bootstrap_result()`](https://morimotoosamu.github.io/lingamr/reference/as_bootstrap_result.md)
は `n_repeats` 次元を（中央値または平均で）集約し、通常の
`BootstrapResult` に変換
するため、既存のブートストラップ照会関数をそのまま適用できる。

``` r

bs_na <- as_bootstrap_result(bwi, aggregate = "median")
get_probabilities(bs_na)
#>    x0 x1 x2 x3 x4 x5
#> x0  0  0  0  1  0  0
#> x1  1  0  1  0  0  0
#> x2  0  0  0  1  0  0
#> x3  0  0  0  0  0  0
#> x4  1  0  1  0  0  0
#> x5  1  0  0  0  0  0
```

[`bootstrap_with_imputation()`](https://morimotoosamu.github.io/lingamr/reference/bootstrap_with_imputation.md)
は総効果を計算しないため、この `BootstrapResult` では
[`get_total_causal_effects()`](https://morimotoosamu.github.io/lingamr/reference/get_total_causal_effects.md)
は利用できない。

## 関連記事

- [手法選択ガイド](https://morimotoosamu.github.io/lingamr/articles/method-selection-ja.md)
  — データに合う手法の選び方
- [Direct LiNGAM
  詳説](https://morimotoosamu.github.io/lingamr/articles/direct-lingam-ja.md)
- [ブートストラップと診断](https://morimotoosamu.github.io/lingamr/articles/bootstrap-diagnostics-ja.md)
