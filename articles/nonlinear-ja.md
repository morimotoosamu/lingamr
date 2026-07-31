# 非線形の手法: RESITとCAM-UV

LiNGAMは因果関係が**線形**であることを仮定する。原因が非線形な関数（飽和、閾値、
2次の用量反応など）を通じて結果に作用する場合、線形の手法は何の警告もなく誤った
方向を選びうる。この記事では `lingamr` の2つの非線形手法を扱う。

- **RESIT**（[`lingam_resit()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit.md)）:
  非線形加法ノイズモデル。潜在交絡変数は*ない*ことを 仮定する。
- **CAM-UV**（[`lingam_camuv()`](https://morimotoosamu.github.io/lingamr/reference/lingam_camuv.md)）:
  加えて**未観測変数**を許容する因果加法モデル。

どちらの手法も、係数ではなく0/1のエッジ指示子を返す。因果関数が非線形である以上、
「$`x_j`$ が $`x_i`$
に及ぼす効果」を1つの数値では要約できないため、総効果の推定は
意図的に提供されていない。

また、どちらも（GAM回帰のために）Suggestsパッケージの `mgcv`
を必要とする。

``` r

library(lingamr)
```

## 非線形データで線形LiNGAMが失敗する例

[`generate_resit_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_resit_sample.md)
は、因果関数がすべて非線形（$`x_1 = 3x_0^2 + e`$、
$`x_2 = 2\tanh(x_1) + 0.8x_0^3 + e`$、$`x_3 = 1.5\sin(x_2) + e`$）の4変数連鎖を生成する。

``` r

nonlin <- generate_resit_sample(n = 300, seed = 1)
head(nonlin$data)
#>           x0          x1           x2           x3
#> 1 -0.4689827  0.83354648  1.596482936  1.831404164
#> 2 -0.2557522 -0.20891458  0.003539709  0.272152307
#> 3  0.1457067  0.05628747 -0.237588016 -0.580258311
#> 4  0.8164156  1.96115504  2.607512659  0.451737357
#> 5 -0.5966361  0.94314057  1.779071173  1.193345462
#> 6  0.7967794  2.39567131  2.846533624 -0.001835138

# True 0/1 adjacency matrix (row = to, column = from)
nonlin$adjacency_matrix
#>    x0 x1 x2 x3
#> x0  0  0  0  0
#> x1  1  0  0  0
#> x2  1  1  0  0
#> x3  0  0  1  0
```

このデータに**線形**のDirect
LiNGAMを適用すると結果は悪くなる。推定された順序は 真の順序
`x0, x1, x2, x3` と一致するとは限らない。

``` r

linear_fit <- lingam_direct(nonlin$data)
colnames(nonlin$data)[linear_fit$causal_order]
#> [1] "x3" "x2" "x1" "x0"
```

これはDirect
LiNGAMの欠陥ではない。単に線形性の仮定がここでは成り立たないという
ことである。非線形の手法が必要になる。

## RESIT：非線形加法ノイズモデル

**RESIT**（Regression with Subsequent Independence Test; Peters et
al. 2014）は、 非線形加法ノイズモデル

``` math
x_i = f_i(\mathrm{parents}(x_i)) + e_i
```

を仮定し、2つのフェーズで構造を復元する。

1.  **順序探索**: 最もsinkらしい変数 –
    残りすべての変数に対する回帰残差が、それら
    と最も独立に近い（HSIC統計量が最小の）変数 – を繰り返し切り離す。
2.  **枝刈り**:
    各候補親について、他の親で回帰した残差がすでに親集合と独立か （HSIC
    p値が `alpha` を上回るか）を検定し、独立なら親から外す。

``` r

resit_result <- lingam_resit(nonlin$data)
print(resit_result)
#> RESIT Result
#>   Variables : 4
#>   Regressor : gam
#>   Causal order: x0 -> x1 -> x2 -> x3
#> 
#> Adjacency matrix (row = to, col = from):
#>   (entries are 0/1 edge indicators, not coefficients)
#>    x0 x1 x2 x3
#> x0  0  0  0  0
#> x1  1  0  0  0
#> x2  1  1  0  0
#> x3  0  0  1  0
```

推定された因果順序と0/1隣接行列は、真の非線形構造を復元している。

``` r

colnames(nonlin$data)[resit_result$causal_order]
#> [1] "x0" "x1" "x2" "x3"
resit_result$adjacency_matrix
#>    x0 x1 x2 x3
#> x0  0  0  0  0
#> x1  1  0  0  0
#> x2  1  1  0  0
#> x3  0  0  1  0
```

エントリはエッジ指示子なので、[`plot_adjacency()`](https://morimotoosamu.github.io/lingamr/reference/plot_adjacency.md)
は重みなしのグラフを描く。

``` r

resit_result$adjacency_matrix |>
  plot_adjacency(
    labels = colnames(nonlin$data),
    title  = "RESIT estimate (edges are 0/1, not coefficients)"
  )
```

### alphaの選択

`alpha`（デフォルト0.01）の向きに注意すること。親が外されるのはHSIC
p値が `alpha` を*上回った*ときなので、`alpha`
を**大きく**するほど独立性の判定が厳しくなり、
結果として**エッジは多く残る**。

### カスタム回帰器

内部の回帰にはデフォルトで平滑化スプラインGAM（`mgcv`）が使われる。当てはめ値を
返す関数 `function(X, y)`
を渡せば、任意の非線形回帰器に差し替えられる。たとえば
ランダムフォレストを使うには次のようにする。

``` r

rf_regressor <- function(X, y) {
  fit <- randomForest::randomForest(X, y)
  as.numeric(predict(fit, X))
}
lingam_resit(nonlin$data, regressor = rf_regressor)
```

### ブートストラップ

[`lingam_resit_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit_bootstrap.md)
は他のブートストラップ関数と同様に機能する。集約された
「確率」は各0/1エッジの検出頻度である。RESITは*リサンプルごとに*
$`O(p^2)`$ 回の GAM当てはめとHSIC検定（各
$`O(n^2)`$）を実行するため、`n_sampling` は控えめにする。

``` r

resit_bs <- lingam_resit_bootstrap(nonlin$data,
  n_sampling = 20L, seed = 1, verbose = FALSE
)
get_probabilities(resit_bs)
#>      [,1] [,2] [,3] [,4]
#> [1,] 0.00 0.00  0.0 0.00
#> [2,] 1.00 0.00  0.1 0.05
#> [3,] 1.00 0.90  0.0 0.00
#> [4,] 0.05 0.05  1.0 0.00
```

## CAM-UV：未観測変数を持つ非線形モデル

RESITは、関連する変数がすべて観測されていることを仮定する。**CAM-UV**
（Maeda and Shimizu
2021）はその仮定を外したものである。変数の部分集合を走査して
残差の独立性を検定することで各変数の直接の親を特定し、特定されたどのエッジでも
説明できない依存が残る変数ペアを、未観測変数を介して結ばれたペアとして報告する。

- **UCP**（unobserved causal path）: 未観測の中間変数を通る有向パス。
- **UBP**（unobserved backdoor path）:
  未観測の共通祖先（潜在交絡変数）。

こうしたペアは、誤っているかもしれない向きを無理に与えられるのではなく、隣接行列
上で `NA` としてマークされる。

[`generate_camuv_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_camuv_sample.md)
はPython版チュートリアルの構造に従う。6つの観測変数の
うち、`u1`（潜在）が `x3` と `x4` を交絡し（UBP）、`u2`（潜在）が `x2`
から `x5` へのパスを媒介する（UCP）。

``` r

camuv_dat <- generate_camuv_sample(n = 500, seed = 1)
head(camuv_dat$data)
#>           x0          x1         x2         x3        x4        x5
#> 1 -0.8276528  0.13274564  0.1053562 2.06976237 2.2118223 0.7967904
#> 2 -0.4513472  0.48681681  0.6321648 1.68111156 1.7569962 0.4177104
#> 3  0.2571408  0.94440653 -0.3991332 0.04983659 0.6192421 1.3267213
#> 4  1.4407966  2.52524438  1.5559126 0.33134200 1.9471524 0.3952331
#> 5 -1.0529335 -0.03798973 -1.3050977 2.38726402 0.9939968 1.8703122
#> 6  1.4061430  2.96254650 -1.5761295 1.77541198 1.4811392 3.7850018

# The UBP (x3-x4) and UCP (x2-x5) entries of the true matrix are NA
camuv_dat$adjacency_matrix
#>    x0 x1 x2 x3 x4 x5
#> x0  0  0  0  0  0  0
#> x1  1  0  0  0  0  0
#> x2  0  0  0  0  0 NA
#> x3  1  0  0  0 NA  0
#> x4  0  0  1 NA  0  0
#> x5  0  0 NA  0  0  0
camuv_dat$confounded_pairs
#>      var1 var2
#> [1,]    3    6
#> [2,]    4    5
```

``` r

camuv_result <- lingam_camuv(camuv_dat$data)
print(camuv_result)
#> CAM-UV Result
#>   Variables : 6
#>   Regressor : gam
#> 
#> Parent sets:
#>   P(x0) = {}
#>   P(x1) = {x0}
#>   P(x2) = {}
#>   P(x3) = {x0}
#>   P(x4) = {x2}
#>   P(x5) = {}
#> 
#> Pairs with an unobserved causal/backdoor path (UCP/UBP):
#>   x2 -- x5
#>   x3 -- x4
#> 
#> Adjacency matrix (row = to, col = from):
#>   (entries are 0/1 edge indicators, not coefficients;
#>    NA = pair connected through an unobserved variable)
#>    x0 x1 x2 x3 x4 x5
#> x0  0  0  0  0  0  0
#> x1  1  0  0  0  0  0
#> x2  0  0  0  0  0 NA
#> x3  1  0  0  0 NA  0
#> x4  0  0  1 NA  0  0
#> x5  0  0 NA  0  0  0
```

[`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)
と同様、因果順序は存在しない。結果は各変数の直接の親と、フラグ
付けされたペアを保持する。

``` r

camuv_result$parents_list
#> $x0
#> integer(0)
#> 
#> $x1
#> [1] 1
#> 
#> $x2
#> integer(0)
#> 
#> $x3
#> [1] 1
#> 
#> $x4
#> [1] 3
#> 
#> $x5
#> integer(0)
camuv_result$confounded_pairs
#>      var1 var2
#> [1,]    3    6
#> [2,]    4    5
camuv_result$adjacency_matrix
#>    x0 x1 x2 x3 x4 x5
#> x0  0  0  0  0  0  0
#> x1  1  0  0  0  0  0
#> x2  0  0  0  0  0 NA
#> x3  1  0  0  0 NA  0
#> x4  0  0  1 NA  0  0
#> x5  0  0 NA  0  0  0
```

### チューニング

- `num_explanatory_vals`（デフォルト2）は、親を探すために走査する変数部分集合の
  サイズの上限。大きくすると統計的検出力は上がるが、コストは組合せ的に増加する。
- `independence = "fcorr"` は独立性の尺度をHSIC検定からF相関（閾値
  `ind_corr`）に 切り替える。`num_explanatory_vals = 2`
  の場合のみサポートされる。
- 事前知識はペア
  `c(i, j)`（「変数iは変数jの原因になりえない」の意、1始まりの
  インデックス）として与える。

Python実装と同様、ブートストラップ版は存在しない。

## RESITとCAM-UVのどちらを使うか

| 観点 | RESIT | CAM-UV |
|----|----|----|
| 潜在交絡変数 | 許容しない | 検出して報告する（UCP/UBPペア） |
| 出力 | 因果順序 + 0/1隣接行列 | 親リスト + `NA` ペアを含む隣接行列 |
| Bootstrap | [`lingam_resit_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit_bootstrap.md) | なし |
| コスト | $`O(p^2)`$ 回のGAM当てはめ + HSIC検定 | `num_explanatory_vals` について組合せ的 |

主要な共通原因がすべて観測されていると主張できるならRESITから始める。潜在交絡が
ありそうならCAM-UVを使う。どちらもHSIC検定ごとに $`n \times n`$
のグラム行列を構築 するため、$`n`$
が数千のオーダーでは推奨されない。必要なら先にサブサンプリングする。

関係が線形で潜在交絡変数が疑われる場合は、
[潜在交絡変数の記事](https://morimotoosamu.github.io/lingamr/articles/latent-confounders-ja.md)の線形手法（ParceLiNGAM、RCD）の
方が計算が軽く、係数の推定値も得られる。

## 関連記事

- [手法選択ガイド](https://morimotoosamu.github.io/lingamr/articles/method-selection-ja.md)
  — データに合う手法の選び方
- [潜在交絡変数（ParceLiNGAMとRCD）](https://morimotoosamu.github.io/lingamr/articles/latent-confounders-ja.md)
  — 線形版の対応手法
- [ブートストラップと診断](https://morimotoosamu.github.io/lingamr/articles/bootstrap-diagnostics-ja.md)
