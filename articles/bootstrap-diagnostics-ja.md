# ブートストラップ安定性とモデル診断

因果探索の推定結果は、その背後にある仮定が成り立つ範囲でしか信頼できない。この
記事では、推定された構造を評価するために `lingamr`
が提供するツール群を扱う。

- **仮定の確認**: 残差の独立性と非ガウス性
  （[`get_error_independence_p_values()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values.md)、[`test_residual_normality()`](https://morimotoosamu.github.io/lingamr/reference/test_residual_normality.md)、
  [`summary_lingam()`](https://morimotoosamu.github.io/lingamr/reference/summary_lingam.md)）。
- **ブートストラップ安定性**:
  推定されたエッジや因果順序はどれだけ再現されるか
  （[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md)
  とその照会関数群）。
- **モデル適合度**:
  推定されたグラフのSEM適合度指標（[`evaluate_model_fit()`](https://morimotoosamu.github.io/lingamr/reference/evaluate_model_fit.md)）。
- **broom連携**: 下流の分析のための
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) /
  [`glance()`](https://generics.r-lib.org/reference/glance.html)。

例では標準の6変数サンプルデータに対するDirect
LiNGAMを使うが、同じワークフローは
他の推定器にも当てはまる（各手法のブートストラップ版はそれぞれの記事で解説
している）。

``` r

library(lingamr)

x1k <- generate_lingam_sample_6(n = 1000)
model <- lingam_direct(x1k$data)
```

## 誤差変数間の独立性

LiNGAMは、残差が独立であると仮定する。[`get_error_independence_p_values()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values.md)
は残差間の 独立性検定のp値を返す。

``` r

p_vals <- x1k$data |>
  get_error_independence_p_values(model)

round(p_vals, 3)
#>       x0    x1    x2    x3    x4    x5
#> x0    NA 0.988 0.214 0.976 0.876 0.952
#> x1 0.988    NA 0.986 0.991 0.328 0.882
#> x2 0.214 0.986    NA 0.919 0.051 0.124
#> x3 0.976 0.991 0.919    NA 0.934 0.978
#> x4 0.876 0.328 0.051 0.934    NA 0.650
#> x5 0.952 0.882 0.124 0.978 0.650    NA
```

p値が小さい（残差間に依存がある）場合、潜在交絡変数や構造の誤指定が示唆される
–
[潜在交絡変数の記事](https://morimotoosamu.github.io/lingamr/articles/latent-confounders-ja.md)の手法を検討すること。

## 残差の正規性検定

残差の正規性を検定する。LiNGAMは非ガウス性を仮定するため、正規性が**棄却される**
（p値が小さい）ことは、モデルの仮定と整合的である。

``` r

# Shapiro-Wilk (default)
x1k$data |>
  test_residual_normality(model)
#> === Residual Normality Test ===
#> Method:         shapiro
#> Sample size:    1000
#> Significance:   0.050
#> Non-Gaussian:   6 / 6 variables
#> 
#>  variable statistic   p_value is_non_gauss skewness kurtosis
#>        x0    0.9516 < 2.2e-16         TRUE    0.061   -1.215
#>        x1    0.9521 < 2.2e-16         TRUE    0.026   -1.213
#>        x2    0.9557 < 2.2e-16         TRUE    0.083   -1.170
#>        x3    0.9578  2.25e-16         TRUE    0.025   -1.163
#>        x4    0.9544 < 2.2e-16         TRUE   -0.003   -1.206
#>        x5    0.9536 < 2.2e-16         TRUE   -0.052   -1.206
#> 
#> Interpretation:
#>   is_non_gauss = TRUE  -> rejects normality (supports LiNGAM assumption)
#>   is_non_gauss = FALSE -> cannot reject normality (LiNGAM may not fit)
#> 
#> All residuals are non-Gaussian. LiNGAM assumption is supported.
```

QQプロットでも残差の正規性を確認する。

``` r

x1k$data |>
  plot_residual_qq(model)
```

![](bootstrap-diagnostics-ja_files/figure-html/qqplot-1.png)

正規性が棄却*されない*場合、因果の方向は識別できていない可能性がある –
[Direct
LiNGAM記事](https://morimotoosamu.github.io/lingamr/articles/direct-lingam-ja.md)の非ガウス性の実験を参照。

## モデルサマリー

[`summary_lingam()`](https://morimotoosamu.github.io/lingamr/reference/summary_lingam.md)
は残差独立性検定と正規性検定をまとめて実行し、LiNGAMが依拠する
2つの仮定（残差が互いに独立であること、残差が非ガウスであること）がどの程度成り立
っているかを一目で確認できるようにする。[`get_error_independence_p_values()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values.md)
と
[`test_residual_normality()`](https://morimotoosamu.github.io/lingamr/reference/test_residual_normality.md)
を個別に呼ぶ代わりに、診断結果を一箇所でまとめて見る ことができる。

``` r

x1k$data |>
  summary_lingam(model)
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

## Direct LiNGAMのブートストラップ

ブートストラップ法を用いてモデルの信頼性を評価する。

``` r

bs_model <- x1k$data |>
  lingam_direct_bootstrap(n_sampling = 100L, seed = 42)
#> Bootstrap: 100 iterations, method=adaptive_lasso (sequential)
#>   iteration 1 / 100
#>   iteration 10 / 100
#>   iteration 20 / 100
#>   iteration 30 / 100
#>   iteration 40 / 100
#>   iteration 50 / 100
#>   iteration 60 / 100
#>   iteration 70 / 100
#>   iteration 80 / 100
#>   iteration 90 / 100
#>   iteration 100 / 100
#> Completed in 3.2 seconds.

bs_model
#> BootstrapResult: 100 samplings, 6 features
```

反復回数や変数の数が多い場合は、`parallel = TRUE`
を指定すると複数コアで高速に 実行できる。コア数は `n_cores`
で指定する（未指定時は安全のため2コアに制限される）。

``` r

bs_model <- x1k$data |>
  lingam_direct_bootstrap(
    n_sampling = 100L,
    seed       = 42,
    parallel   = TRUE,
    n_cores    = 4L
  )
```

並列実行ではL’Ecuyerの並列乱数ストリームが使われるため、同じ `seed`
と同じ `n_cores`
であれば結果は再現可能だが、逐次実行（`parallel = FALSE`）の結果とは
数値的に一致しない点に注意する。

### ブートストラップ結果の確認

ブートストラップ結果から、各パスの出現頻度と係数の平均を計算する。

``` r

bs_model |>
  get_causal_direction_counts(labels = names(x1k$data))
#>    from to count proportion mean_effect median_effect  sd_effect    ci_lower
#> 1     1  6   100       1.00  4.01532920    4.01513886 0.01126767  3.99550980
#> 2     1  2    99       0.99  2.98181621    2.97864538 0.02849338  2.92980702
#> 3     1  5    99       0.99  8.00994011    8.00748238 0.02951185  7.95680521
#> 4     3  2    99       0.99  2.00498455    2.00660933 0.01479861  1.97675886
#> 5     3  5    99       0.99 -1.00529230   -1.00485827 0.01523485 -1.03801290
#> 6     4  1    99       0.99  3.03521019    3.03586526 0.03001961  2.97855949
#> 7     4  3    99       0.99  5.99644109    5.99745219 0.03186571  5.94050363
#> 8     2  1     1       0.01  0.05304916    0.05304916 0.00000000  0.05304916
#> 9     2  3     1       0.01  0.40196452    0.40196452 0.00000000  0.40196452
#> 10    2  5     1       0.01  0.90679690    0.90679690 0.00000000  0.90679690
#> 11    3  4     1       0.01  0.16166764    0.16166764 0.00000000  0.16166764
#> 12    5  1     1       0.01  0.10453910    0.10453910 0.00000000  0.10453910
#> 13    5  3     1       0.01 -0.13636255   -0.13636255 0.00000000 -0.13636255
#>       ci_upper from_name to_name
#> 1   4.03698551        x0      x5
#> 2   3.03860193        x0      x1
#> 3   8.07414013        x0      x4
#> 4   2.03193816        x2      x1
#> 5  -0.97488336        x2      x4
#> 6   3.09304642        x3      x0
#> 7   6.06134091        x3      x2
#> 8   0.05304916        x1      x0
#> 9   0.40196452        x1      x2
#> 10  0.90679690        x1      x4
#> 11  0.16166764        x2      x3
#> 12  0.10453910        x4      x0
#> 13 -0.13636255        x4      x2
```

### 平均因果効果の隣接行列

ブートストラップ結果から隣接行列を構築する。

``` r

bs_adjacency_matrix <- bs_model |>
  get_adjacency_matrix_summary(stat = "median")

bs_adjacency_matrix |>
  round(3)
#>       [,1]  [,2]   [,3]  [,4]   [,5] [,6]
#> [1,] 0.000 0.053  0.000 3.036  0.105    0
#> [2,] 2.979 0.000  2.007 0.000  0.000    0
#> [3,] 0.000 0.402  0.000 5.997 -0.136    0
#> [4,] 0.000 0.000  0.162 0.000  0.000    0
#> [5,] 8.007 0.907 -1.005 0.000  0.000    0
#> [6,] 4.015 0.000  0.000 0.000  0.000    0
```

推定した隣接行列を可視化する。

``` r

bs_adjacency_matrix |>
  plot_adjacency(
    labels    = colnames(x1k$data),
    title     = "Estimated (with Bootstrap)",
    rankdir   = "TB",
    shape     = "circle",
    fillcolor = "lightgreen"
  )
```

### パス出現頻度の行列

各パスの出現頻度行列を計算する。

``` r

bs_model |>
  get_probabilities()
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,] 0.00 0.01 0.00 0.99 0.01    0
#> [2,] 0.99 0.00 0.99 0.00 0.00    0
#> [3,] 0.00 0.01 0.00 0.99 0.01    0
#> [4,] 0.00 0.00 0.01 0.00 0.00    0
#> [5,] 0.99 0.01 0.99 0.00 0.00    0
#> [6,] 1.00 0.00 0.00 0.00 0.00    0
```

### 平均総効果

各パスの平均総効果を計算する。

``` r

bs_model |>
  get_total_causal_effects()
#>    from to      effect probability
#> 1     1  6  4.01520158        1.00
#> 2     1  2  2.87431611        0.99
#> 3     1  5  7.90813117        0.99
#> 4     3  2  1.95874622        0.99
#> 5     3  5 -1.06193484        0.99
#> 6     4  1  3.03586526        0.99
#> 7     4  2 21.07027271        0.99
#> 8     4  3  5.99805118        0.99
#> 9     4  5 18.28272145        0.99
#> 10    4  6 12.18719857        0.99
#> 11    3  6 -0.24574320        0.04
#> 12    2  1  0.14794503        0.01
#> 13    2  3  0.27850920        0.01
#> 14    2  4  0.04611007        0.01
#> 15    2  5  0.90679690        0.01
#> 16    2  6  0.59359217        0.01
#> 17    3  4  0.16192779        0.01
#> 18    5  1  0.10498716        0.01
#> 19    5  3 -0.13625059        0.01
#> 20    5  6  0.42156703        0.01
#> 21    6  2  0.24518629        0.01
```

ブートストラップ結果を因果グラフに変換する。デフォルトでは、サンプルの50%以上で
出現したパスのみが表示される。

``` r

bs_model |>
  plot_bootstrap_probabilities()
```

### 2変数間のパス

[`get_paths()`](https://morimotoosamu.github.io/lingamr/reference/get_paths.md)
は2変数間の総因果効果を、それを媒介する個々のパスとそのブートスト
ラップ確率に分解する。[Direct
LiNGAM記事](https://morimotoosamu.github.io/lingamr/articles/direct-lingam-ja.md)で見たとおり、
[`generate_lingam_sample_6()`](https://morimotoosamu.github.io/lingamr/reference/generate_lingam_sample_6.md)
の真の構造では、x3からx1へのパスが2つ存在する:
`x3 -> x0 -> x1`（間接効果9.0）と
`x3 -> x2 -> x1`（間接効果12.0）。インデックス
は1始まりなので、x3（4列目）からx1（2列目）は次のように指定する。

``` r

bs_model |>
  get_paths(4, 2)
#>      path    effect probability
#> 1 4, 1, 2  9.041187        0.99
#> 2 4, 3, 2 12.020020        0.99
```

いずれのパスも100サンプル中99サンプルで検出され、効果の中央値は真の値（9.0と
12.0）に近い。

### 繰り返し出現するDAG構造の頻度

[`get_directed_acyclic_graph_counts()`](https://morimotoosamu.github.io/lingamr/reference/get_directed_acyclic_graph_counts.md)
は各ブートストラップサンプルで推定された
グラフ全体に注目し、個々のエッジではなく異なるDAGがどれだけ繰り返し出現するかを
数える。`n_dags` で頻度上位何件を返すかを絞り込める。

``` r

dag_counts <- bs_model |>
  get_directed_acyclic_graph_counts(n_dags = 3)

dag_counts$count
#> [1] 99  1

dag_counts$dag[[1]]
#>   from to
#> 1    1  2
#> 2    1  5
#> 3    1  6
#> 4    3  2
#> 5    3  5
#> 6    4  1
#> 7    4  3
```

最も頻度の高いDAG（100サンプル中99サンプル）は、[`generate_lingam_sample_6()`](https://morimotoosamu.github.io/lingamr/reference/generate_lingam_sample_6.md)
の 真のエッジ集合と完全に一致する。

### 因果順序の安定性

[`get_causal_order_stability()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_order_stability.md)
は各ブートストラップサンプルで推定された因果順序を
集約し、その順序がどれだけ安定しているかを定量化する。各変数の順位分布、変数ペアの
先行確率（`P[i, j]` =
変数iが変数jより上流にあったサンプルの割合）、および全体の
安定性スコア（0 = ランダム、1 = 全サンプルで完全に一致）を返す。

``` r

bs_model |>
  get_causal_order_stability(labels = names(x1k$data))
#> === Causal Order Stability ===
#> Bootstrap samples:       100
#> Overall stability score: 0.736  (0 = random, 1 = fully stable)
#> 
#> Rank summary (sorted by mean rank; 1 = most upstream):
#>  variable mean_rank sd_rank median_rank mode_rank
#>        x3      1.05    0.50           1         1
#>        x0      2.62    0.51           3         3
#>        x2      2.75    0.95           2         2
#>        x5      4.41    1.23           4         3
#>        x4      4.92    0.77           5         5
#>        x1      5.25    0.88           5         6
#> 
#> Precedence probability P[i, j] = P(variable i precedes j):
#>      x0   x1   x2   x3   x4   x5
#> x0 0.00 0.99 0.39 0.01 0.99 1.00
#> x1 0.01 0.00 0.01 0.01 0.38 0.34
#> x2 0.61 0.99 0.00 0.01 0.99 0.65
#> x3 0.99 0.99 0.99 0.00 0.99 0.99
#> x4 0.01 0.62 0.01 0.01 0.00 0.43
#> x5 0.00 0.66 0.35 0.01 0.57 0.00
```

> **注意:**
> ブートストラップの安定性は正しさを保証しない。モデルの仮定が破られて
> いる場合、誤った構造が*安定して*再現されることがある – [Direct
> LiNGAM記事](https://morimotoosamu.github.io/lingamr/articles/direct-lingam-ja.md)の測定誤差パラドックスに実例がある。

## モデル適合度の評価

[`evaluate_model_fit()`](https://morimotoosamu.github.io/lingamr/reference/evaluate_model_fit.md)
は、推定された隣接行列を構造方程式モデル（SEM）として扱い、 `lavaan`
パッケージ（オプション依存。`install.packages("lavaan")`
でインストール）
経由で標準的なSEM適合度指標（CFI、RMSEA、AIC/BICなど）を報告する。これは、推定
された因果グラフが、どのように推定されたかとは独立に、データとどれだけ整合的かを
判断する際に有用である。

``` r

sample6 <- generate_lingam_sample_6()
fit_result <- lingam_direct(sample6$data, reg_method = "ols")

# fit measures for the estimated graph
evaluate_model_fit(fit_result, sample6$data)
#>   DoF DoF Baseline chi2 chi2 p-value chi2 Baseline CFI GFI AGFI NFI TLI RMSEA
#> 1   0           15    0           NA       23023.7   1   1   NA   1   1     0
#>        AIC      BIC    LogLik
#> 1 1860.598 1958.753 -910.2991
```

すべてのエッジの方向を逆転させると誤指定のモデルになり、その適合度指標は明らかに
悪化する（CFIは低下し、RMSEAは上昇する）。

``` r

reversed_adjacency <- t(fit_result$adjacency_matrix)
evaluate_model_fit(reversed_adjacency, sample6$data)
#>   DoF DoF Baseline         chi2 chi2 p-value chi2 Baseline CFI GFI AGFI NFI TLI
#> 1   0           15 2.664535e-11           NA       23023.7   1   1   NA   1   1
#>   RMSEA       AIC       BIC   LogLik
#> 1     0 -4264.864 -4166.708 2152.432
```

## broomとの連携（tidy / glance）

推定結果は、`broom` 互換の
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) /
[`glance()`](https://generics.r-lib.org/reference/glance.html)
を使ってdata.frameに変換でき、 `ggplot2` や `dplyr`
との連携が容易になる。[`tidy()`](https://generics.r-lib.org/reference/tidy.html)
はエッジリスト（`from`、`to`、
`estimate`）を返し、[`glance()`](https://generics.r-lib.org/reference/glance.html)
はモデル全体の1行サマリーを返す。[`tidy()`](https://generics.r-lib.org/reference/tidy.html)
はブート
ストラップ結果に対しても機能し、その場合は各方向の出現頻度などを返す。

``` r

# Convert the estimated adjacency matrix to an edge list
tidy(model)
#>   from to  estimate
#> 1   x0 x1  2.987705
#> 2   x0 x4  8.000096
#> 3   x0 x5  4.014962
#> 4   x2 x1  2.001708
#> 5   x2 x4 -1.000306
#> 6   x3 x0  3.032952
#> 7   x3 x2  5.992677

# One-row summary of the whole model
glance(model)
#>   n_variables n_edges                     causal_order
#> 1           6       7 x3 -> x2 -> x0 -> x4 -> x5 -> x1

# Direction-wise summary of the bootstrap results (variable names via labels)
tidy(bs_model, labels = names(x1k$data))
#>    from to count proportion mean_effect median_effect  sd_effect    ci_lower
#> 1     1  6   100       1.00  4.01532920    4.01513886 0.01126767  3.99550980
#> 2     1  2    99       0.99  2.98181621    2.97864538 0.02849338  2.92980702
#> 3     1  5    99       0.99  8.00994011    8.00748238 0.02951185  7.95680521
#> 4     3  2    99       0.99  2.00498455    2.00660933 0.01479861  1.97675886
#> 5     3  5    99       0.99 -1.00529230   -1.00485827 0.01523485 -1.03801290
#> 6     4  1    99       0.99  3.03521019    3.03586526 0.03001961  2.97855949
#> 7     4  3    99       0.99  5.99644109    5.99745219 0.03186571  5.94050363
#> 8     2  1     1       0.01  0.05304916    0.05304916 0.00000000  0.05304916
#> 9     2  3     1       0.01  0.40196452    0.40196452 0.00000000  0.40196452
#> 10    2  5     1       0.01  0.90679690    0.90679690 0.00000000  0.90679690
#> 11    3  4     1       0.01  0.16166764    0.16166764 0.00000000  0.16166764
#> 12    5  1     1       0.01  0.10453910    0.10453910 0.00000000  0.10453910
#> 13    5  3     1       0.01 -0.13636255   -0.13636255 0.00000000 -0.13636255
#>       ci_upper from_name to_name
#> 1   4.03698551        x0      x5
#> 2   3.03860193        x0      x1
#> 3   8.07414013        x0      x4
#> 4   2.03193816        x2      x1
#> 5  -0.97488336        x2      x4
#> 6   3.09304642        x3      x0
#> 7   6.06134091        x3      x2
#> 8   0.05304916        x1      x0
#> 9   0.40196452        x1      x2
#> 10  0.90679690        x1      x4
#> 11  0.16166764        x2      x3
#> 12  0.10453910        x4      x0
#> 13 -0.13636255        x4      x2
```

## 関連記事

- [手法選択ガイド](https://morimotoosamu.github.io/lingamr/articles/method-selection-ja.md)
- [Direct LiNGAM
  詳説](https://morimotoosamu.github.io/lingamr/articles/direct-lingam-ja.md)
  — ブートストラップが安定して誤る 測定誤差パラドックスを含む
- [時系列](https://morimotoosamu.github.io/lingamr/articles/time-series-ja.md)、[潜在交絡変数](https://morimotoosamu.github.io/lingamr/articles/latent-confounders-ja.md)、
  [特殊なデータ](https://morimotoosamu.github.io/lingamr/articles/special-data-ja.md)
  — 各手法のブートストラップ版
