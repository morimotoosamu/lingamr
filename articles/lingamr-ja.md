# lingamrによる因果探索

`lingamr` は、LiNGAM 系のアルゴリズムを使って、純粋な観測データから
**因果構造**（どの変数がどの変数の原因であり、どれだけ強く影響するか）を推定する
（[清水研究室](https://www.shimizulab.org/)による Python
[`lingam`](https://github.com/cdt15/lingam) パッケージのR移植版）。

この vignette
では、LiNGAMの中核となる考え方を説明し、最小限のエンドツーエンドの
ワークフローを一通り実行する。各手法の詳細なガイドはパッケージサイトにある
（[次に読むべきもの](#where-to-go-next) を参照）。

``` r

library(lingamr)
```

## 基本的な考え方：なぜ非ガウス性から因果の方向がわかるのか

相関だけでは「xがyの原因」なのか「yがxの原因」なのかを区別できない。どちらの
モデルもまったく同じ共分散行列を生成しうるからである。したがって、2次モーメント
（分散・共分散）に基づく古典的な手法が返せるのは、せいぜい構造の*同値類*までである。

LiNGAM（Linear Non-Gaussian Acyclic Model; Shimizu et
al. 2006）は、誤差項が
**非ガウス**であるという仮定を1つ加えることで方向を解決する。各変数のモデルは

``` math
x_i = \sum_{j:\ \mathrm{parent\ of\ } i} b_{ij}\, x_j + e_i
```

であり、誤差 $`e_i`$
は互いに独立な非ガウス分布に従い、構造はDAG（有向非巡回グラフ）
をなす。これらの仮定の下で、因果構造は観測データから同値類ではなく**一意に識別可能**
になる。

直感は1つの実験で確認できる。真のモデルを
$`y = 1.5x + e`$（誤差は*一様分布*、
すなわち非ガウス）とし、両方向に回帰してみる。

``` r

n <- 1000
x <- runif(n, -1, 1)
y <- 1.5 * x + runif(n, -1, 1)

res_causal <- residuals(lm(y ~ x)) # correct direction:  x -> y
res_anti   <- residuals(lm(x ~ y)) # reverse direction:  y -> x

oldpar <- par(mfrow = c(1, 2))
plot(x, res_causal, main = "Correct: residual of y ~ x", cex = 0.4)
plot(y, res_anti,   main = "Reverse: residual of x ~ y", cex = 0.4)
```

![](lingamr-ja_files/figure-html/nongauss_intuition-1.png)

``` r

par(oldpar)
```

正しい方向では、残差は説明変数と**独立**である（左パネル：のっぺりした帯）。逆方向
ではそうならない（右パネル：残差の散らばりが $`y`$ に依存する）。Direct
LiNGAMはこの
非対称性をアルゴリズムに変換したものである。残差が最も独立になる変数こそが最上流
（外生）の変数であり、それを取り除き、回帰で影響を除去し、同じ手続きを繰り返す。

誤差がガウス分布であれば、両パネルは同じに見えてしまう。これが非ガウス性が不可欠な
理由であり、*基本の* LiNGAMモデルが次の4つの仮定に立脚する理由である。

1.  **線形**な関係
2.  **非巡回**なグラフ（DAG）
3.  **非ガウス**で互いに独立な誤差
4.  **潜在交絡変数がない**（すべての共通原因が観測されている）

これに加えて観測はi.i.d.であることを仮定する。`lingamr`
には、これらの仮定を
それぞれ緩和した推定器も収録されている（[`vignette("method-selection-ja")`](https://morimotoosamu.github.io/lingamr/articles/method-selection-ja.md)
を参照）。

## 最小限のワークフロー

### 推定

[`generate_lingam_sample_6()`](https://morimotoosamu.github.io/lingamr/reference/generate_lingam_sample_6.md)
は既知の6変数LiNGAMモデルからデータを生成するので、
推定結果を真の構造と比較できる。[`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
はDirect LiNGAMを実行する。
デフォルトでは独立性は相互情報量で評価され、係数はadaptive
LASSOで推定される。

``` r

x1k <- generate_lingam_sample_6(n = 1000)

model <- lingam_direct(x1k$data)
model
#> Direct LiNGAM Result
#>   Variables : 6
#>   Causal order: x3 -> x2 -> x0 -> x4 -> x5 -> x1
#> 
#> Adjacency matrix (row = to, col = from):
#>       x0 x1     x2    x3 x4 x5
#> x0 0.000  0  0.000 3.033  0  0
#> x1 2.988  0  2.002 0.000  0  0
#> x2 0.000  0  0.000 5.993  0  0
#> x3 0.000  0  0.000 0.000  0  0
#> x4 8.000  0 -1.000 0.000  0  0
#> x5 4.015  0  0.000 0.000  0  0
```

結果の2つの主要な要素:

``` r

# 推定された因果順序（上流が先頭）
colnames(x1k$data)[model$causal_order]
#> [1] "x3" "x2" "x0" "x4" "x5" "x1"

# 隣接行列: B[i, j] は x_j から x_i への直接効果
round(model$adjacency_matrix, 3)
#>       x0 x1     x2    x3 x4 x5
#> x0 0.000  0  0.000 3.033  0  0
#> x1 2.988  0  2.002 0.000  0  0
#> x2 0.000  0  0.000 5.993  0  0
#> x3 0.000  0  0.000 0.000  0  0
#> x4 8.000  0 -1.000 0.000  0  0
#> x5 4.015  0  0.000 0.000  0  0
```

### 可視化

[`plot_adjacency()`](https://morimotoosamu.github.io/lingamr/reference/plot_adjacency.md)
は因果グラフを描画する。`true_B` に真の構造を渡すと、比較結果が
色分けされる（緑 = 正解、赤 = 誤検出、オレンジ破線 = 見逃し）。

``` r

model$adjacency_matrix |>
  plot_adjacency(
    labels = colnames(x1k$data),
    true_B = x1k$true_adjacency,
    title  = "Estimated vs. true structure"
  )
```

### 介入の効果：総因果効果

総因果効果とは、ある変数を1単位動かしたとき、すべてのパスを通じて別の変数が最終的に
どれだけ動くかであり、介入について考えるにはこちらが必要になる（重回帰係数は別の
問いに答えるものである。詳細は [Direct
LiNGAM記事](https://morimotoosamu.github.io/lingamr/articles/direct-lingam-ja.html)
を参照）。

``` r

total_effects <- estimate_all_total_effects(x1k$data, model)
round(total_effects, 2)
#>      x0 x1    x2    x3 x4 x5
#> x0 0.00  0  0.00  3.03  0  0
#> x1 2.87  0  1.94 21.06  0  0
#> x2 0.00  0  0.00  5.99  0  0
#> x3 0.00  0  0.00  0.00  0  0
#> x4 7.91  0 -1.13 18.28  0  0
#> x5 4.02  0  0.00 12.18  0  0
```

### 仮定の確認

[`summary_lingam()`](https://morimotoosamu.github.io/lingamr/reference/summary_lingam.md)
は2つの主要な診断をまとめて実行する。残差は互いに**独立**で
あるべきで（仮定4）、かつ**非ガウス**であるべきである（仮定3。したがって、ここでは
正規性が*棄却される*ことが良い知らせになる）。

``` r

summary_lingam(x1k$data, model)
#> === Direct LiNGAM Model Summary ===
#> Variables:    6
#> Observations: 1000
#> Edges:        7
#> Causal order: x3 -> x2 -> x0 -> x4 -> x5 -> x1
#> 
#> --- Assumption 1: Independence of residuals ---
#> Method:           spearman
#> Dependent pairs:  0 / 15  (p < 0.050)
#> Min p-value:      0.0510
#> => Residuals appear mutually independent (assumption supported).
#> 
#> --- Assumption 2: Non-Gaussianity of residuals ---
#> Method:           shapiro
#> Non-Gaussian:     6 / 6  (p <= 0.050)
#> => All residuals are non-Gaussian (assumption supported).
```

### 安定性の定量化

ブートストラップはリサンプルしたデータで推定を再実行し、各エッジがどれだけの頻度で
再現されるかを報告する。出現確率の低いエッジを過大解釈してはならない。

``` r

bs <- lingam_direct_bootstrap(x1k$data, n_sampling = 50L, seed = 42)
#> Bootstrap: 50 iterations, method=adaptive_lasso (sequential)
#>   iteration 1 / 50
#>   iteration 10 / 50
#>   iteration 20 / 50
#>   iteration 30 / 50
#>   iteration 40 / 50
#>   iteration 50 / 50
#> Completed in 1.3 seconds.

round(get_probabilities(bs), 2)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,] 0.00 0.02 0.00 0.98 0.02    0
#> [2,] 0.98 0.00 0.98 0.00 0.00    0
#> [3,] 0.00 0.02 0.00 0.98 0.02    0
#> [4,] 0.00 0.00 0.02 0.00 0.00    0
#> [5,] 0.98 0.02 0.98 0.00 0.00    0
#> [6,] 1.00 0.00 0.00 0.00 0.00    0
```

## 次に読むべきもの

時系列・潜在交絡・非線形関係・混合データ・欠測データなど、あなたのデータに合う
推定器の選び方は
[`vignette("method-selection-ja")`](https://morimotoosamu.github.io/lingamr/articles/method-selection-ja.md)
で解説している。

詳細な実行例は[パッケージサイト](https://morimotoosamu.github.io/lingamr/)にある。

| 記事 | 内容 |
|----|----|
| [Direct LiNGAM 詳説](https://morimotoosamu.github.io/lingamr/articles/direct-lingam-ja.html) | 事前知識、回帰手法、非ガウス性の実験、高次元データ、失敗するケース |
| [ブートストラップと診断](https://morimotoosamu.github.io/lingamr/articles/bootstrap-diagnostics-ja.html) | 安定性分析、仮定の確認、SEM適合度、broom連携 |
| [時系列](https://morimotoosamu.github.io/lingamr/articles/time-series-ja.html) | VAR-LiNGAM、VARMA-LiNGAM |
| [潜在交絡変数](https://morimotoosamu.github.io/lingamr/articles/latent-confounders-ja.html) | BottomUpParceLiNGAM、RCD |
| [非線形の手法](https://morimotoosamu.github.io/lingamr/articles/nonlinear-ja.html) | RESIT、CAM-UV |
| [特殊なデータ](https://morimotoosamu.github.io/lingamr/articles/special-data-ja.html) | 混合データ（LiM）、複数グループ、欠測値 |

英語版は
[`vignette("lingamr")`](https://morimotoosamu.github.io/lingamr/articles/lingamr.md)
を参照。サイトの記事にも英語版がある。
