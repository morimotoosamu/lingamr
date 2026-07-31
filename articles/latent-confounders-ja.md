# 潜在交絡変数: ParceLiNGAMとRCD

標準のDirect
LiNGAMは、潜在（未観測の）交絡変数が存在しないことを仮定している。
すなわち、観測変数のうち2つ以上に影響を及ぼす変数は、それ自体も観測されていなければ
ならない。この仮定が成り立たない場合、[`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
は依然として完全な因果順序を
返すが、それは無言のまま行われる。その一部が誤っている可能性があっても、どの部分が
誤っているかは示されない。

この記事では、潜在交絡変数に対応する**線形**の2手法を扱う。

- **BottomUpParceLiNGAM**（[`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)）:
  順序付けできなかった変数を単一の
  *未解決ブロック*として返し、潜在交絡変数がありそうな箇所を示す。
- **RCD**（[`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)）:
  各変数の祖先集合を推定し、潜在交絡変数を共有する
  *特定のペア*を明示する。

因果関係が**非線形**でもありうる場合は、
[非線形手法の記事](https://morimotoosamu.github.io/lingamr/articles/nonlinear-ja.md)のCAM-UVを参照。非線形モデルにおける未観測
変数を扱う手法である。

``` r

library(lingamr)
```

## Bottom-Up ParceLiNGAM

[`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)（BottomUpParceLiNGAM、Tashiro
et al. 2014）は、sink（最も下流）側
から因果順序を探索し、各ステップで候補変数の残差が他の変数と独立かどうかを検定する。
その検定が棄却された時点で探索は停止し、まだ配置できていないすべての変数は、単一の
**未解決ブロック**としてまとめて返される –
これは、それらの変数がおそらく潜在交絡
変数を共有していることを示すシグナルであり、順序に関する（誤っているかもしれない）
推測ではない。

[`generate_parce_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_parce_sample.md)
は、`x6` が `x2` と `x3` の未観測の共通原因である7変数
モデルを生成する。データとして返されるのは `x0`-`x5` のみである。

``` r

# HSIC is O(n^2), so a moderate n keeps this article fast to build
confounded <- generate_parce_sample(n = 500, seed = 1)
head(confounded$data)
#>          x0        x1        x2        x3          x4        x5
#> 1 0.6154746 1.7104554 1.0618261 1.0851944  0.57917376 0.8337563
#> 2 1.5905703 2.4770365 1.4291087 1.4325230  0.56022597 0.8686206
#> 3 1.1007549 2.1817888 1.5289901 1.8037643 -0.03671602 1.4001192
#> 4 1.7744689 2.7106515 2.7714036 2.4797583 -0.10133399 1.3102925
#> 5 0.5433612 0.7244786 0.5217204 0.8755981  0.82504738 1.2597767
#> 6 1.7671488 1.8838085 1.8358794 2.7663075  0.90699176 1.3624485
confounded$confounded_pair
#> [1] 3 4
```

``` r

parce_result <- lingam_parce(confounded$data, reg_method = "ols")
print(parce_result)
#> Bottom-Up ParceLiNGAM Result
#>   Variables : 6
#>   Independence measure: hsic
#>   Causal order: (x2, x3) -> x0 -> x4 -> x5 -> x1
#>   (NA entries in the adjacency matrix = unresolved order / suspected latent confounding)
#> 
#> Adjacency matrix (row = to, col = from):
#>       x0 x1     x2     x3    x4     x5
#> x0 0.000  0 -0.010  0.516 0.000  0.000
#> x1 0.479  0  0.447  0.060 0.025 -0.049
#> x2 0.000  0  0.000     NA 0.000  0.000
#> x3 0.000  0     NA  0.000 0.000  0.000
#> x4 0.497  0 -0.490 -0.001 0.000  0.000
#> x5 0.436  0  0.068  0.023 0.050  0.000
```

因果順序の最初の要素は未解決ブロックであり、括弧付きで表示される。ここでは
`x2` と `x3` が正しく含まれている。隣接行列の対応するエントリは `NA`
になり、 残りの完全に解決された変数間のエッジは通常通り推定される。

``` r

parce_result$causal_order[[1]]
#> [1] 3 4
parce_result$adjacency_matrix[confounded$confounded_pair, confounded$confounded_pair]
#>    x2 x3
#> x2  0 NA
#> x3 NA  0
```

交絡された変数の真の親は特定できないため、[`estimate_total_effect_parce()`](https://morimotoosamu.github.io/lingamr/reference/estimate_total_effect_parce.md)
は
未解決ブロック内の変数を*起点*とする総効果を求められた場合には警告を発して
`NA` を 返すが、十分に識別可能なペアについては通常通り推定値を計算する。

``` r

# from a confounded variable: warns and returns NA
estimate_total_effect_parce(confounded$data, parce_result,
  from_index = confounded$confounded_pair[1], to_index = "x1"
)
#> Warning in estimate_total_effect_parce(confounded$data, parce_result,
#> from_index = confounded$confounded_pair[1], : x2 is part of an unresolved
#> causal order (suspected latent confounding); total effect cannot be estimated.
#> [1] NA

# a well-identified pair: a normal numeric estimate
estimate_total_effect_parce(confounded$data, parce_result,
  from_index = "x0", to_index = "x5"
)
#> [1] 0.5121874
```

[`lingam_parce_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce_bootstrap.md)
は、[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md)
と同様のスタイルで ブートストラップの安定性推定を提供する。集約時には
`NA`（未解決）のエッジは存在
しないものとして扱われるため、[`get_probabilities()`](https://morimotoosamu.github.io/lingamr/reference/get_probabilities.md)
などの `BootstrapResult` 照会 関数は通常通り機能する。ただし
[`get_causal_order_stability()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_order_stability.md)
だけは例外で、
ParceLiNGAMのブロック化された因果順序はその固定長フォーマットに合わないためである。

``` r

parce_bs <- lingam_parce_bootstrap(confounded$data,
  n_sampling = 10L, reg_method = "ols", seed = 1, verbose = FALSE
)
get_probabilities(parce_bs)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]  0.0  0.2  0.5  0.4  0.0  0.0
#> [2,]  0.5  0.0  0.5  0.4  0.2  0.3
#> [3,]  0.0  0.0  0.0  0.0  0.0  0.0
#> [4,]  0.1  0.2  0.2  0.0  0.1  0.0
#> [5,]  0.7  0.5  0.7  0.6  0.0  0.3
#> [6,]  0.6  0.4  0.6  0.6  0.4  0.0
```

## RCD（Repetitive Causal Discovery）

[`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)（Repetitive
Causal Discovery; Maeda and Shimizu 2020）は、
[`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)
と同じ潜在交絡変数の問題を扱うが、異なる角度からアプローチする。
検定が棄却された時点で**未解決ブロック**として諦めて因果順序を探索するのではなく、
RCDは各変数の**祖先集合**を直接推定し、そのうえで親を持たない個々のペアについて
潜在交絡変数を共有していないかを検定する。この結果、RCDの出力はブロック単位
（順序付けできなかった変数の集合）ではなく、ペア単位（どの特定のペアが交絡されて
いるか）となる。

[`generate_rcd_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_rcd_sample.md)
は、`x6` が `x2` と `x4` の未観測の共通原因である7変数
モデルを生成する。データとして返されるのは `x0`-`x5` のみである。

``` r

# HSIC is O(n^2), so a moderate n keeps this article fast to build
rcd_confounded <- generate_rcd_sample(n = 300, seed = 1)
head(rcd_confounded$data)
#>            x0           x1         x2           x3         x4            x5
#> 1 -0.96839700 -0.023398020 -0.7514703 -0.473147026 -0.9674898 -0.0307310344
#> 2  1.31480945  0.424388510  1.0389688  0.001304296  1.4111680  0.0007741684
#> 3 -0.85973904 -0.025330472 -1.1731736 -0.034157640 -1.3389435 -0.0729373397
#> 4 -0.76463976  0.324413350 -0.7116352  0.078719765 -0.3333176  0.5074829194
#> 5  0.08152383  0.002364107 -0.2546673  0.036715357 -0.2040856  0.0044720536
#> 6 -0.30389184 -0.225029302 -0.5032249 -0.172267536  0.5221642 -0.0690391705
rcd_confounded$confounded_pair
#> [1] 3 5
```

``` r

rcd_result <- lingam_rcd(rcd_confounded$data)
print(rcd_result)
#> RCD Result
#>   Variables : 6
#> 
#> Ancestor sets:
#>   M(x0) = {x1, x3, x5}
#>   M(x1) = {x5}
#>   M(x2) = {x0, x1, x3, x5}
#>   M(x3) = {x5}
#>   M(x4) = {x0, x1, x3, x5}
#>   M(x5) = {}
#> 
#>   (NA entries in the adjacency matrix = suspected shared latent confounder)
#> 
#> Adjacency matrix (row = to, col = from):
#>       x0    x1 x2    x3 x4    x5
#> x0 0.000 1.116  0 0.989  0 0.000
#> x1 0.000 0.000  0 0.000  0 0.588
#> x2 0.810 0.000  0 0.000 NA 0.000
#> x3 0.000 0.000  0 0.000  0 0.449
#> x4 1.015 0.000 NA 0.000  0 0.000
#> x5 0.000 0.000  0 0.000  0 0.000
```

`ancestors_list`
は各変数の推定された祖先を与える（因果順序ではない）。交絡された
ペアの隣接行列エントリは `NA` になる。

``` r

rcd_result$ancestors_list
#> $x0
#> [1] 2 4 6
#> 
#> $x1
#> [1] 6
#> 
#> $x2
#> [1] 1 2 4 6
#> 
#> $x3
#> [1] 6
#> 
#> $x4
#> [1] 1 2 4 6
#> 
#> $x5
#> integer(0)
rcd_result$adjacency_matrix[rcd_confounded$confounded_pair, rcd_confounded$confounded_pair]
#>    x2 x4
#> x2  0 NA
#> x4 NA  0
```

ParceLiNGAMと同様、[`estimate_total_effect_rcd()`](https://morimotoosamu.github.io/lingamr/reference/estimate_total_effect_rcd.md)
は交絡された変数を起点とする総効果 を求められた場合、警告を発して `NA`
を返す。

``` r

# from a confounded variable: warns and returns NA
estimate_total_effect_rcd(rcd_confounded$data, rcd_result,
  from_index = rcd_confounded$confounded_pair[1], to_index = rcd_confounded$confounded_pair[2]
)
#> Warning in estimate_total_effect_rcd(rcd_confounded$data, rcd_result,
#> from_index = rcd_confounded$confounded_pair[1], : x2 is part of a suspected
#> latent confounder pair; total effect cannot be estimated.
#> [1] NA

# a well-identified pair: a normal numeric estimate
estimate_total_effect_rcd(rcd_confounded$data, rcd_result,
  from_index = "x5", to_index = "x0"
)
#>      x5 
#> 1.05674
```

[`lingam_rcd_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd_bootstrap.md)
はブートストラップの安定性推定を提供し、通常の `BootstrapResult`
照会関数がそのまま使える。

## ParceLiNGAMとRCDのどちらを使うか

| 観点 | ParceLiNGAM | RCD |
|----|----|----|
| 出力の粒度 | ブロック単位: 順序付けできなかった変数の集合 | ペア単位: どの特定のペアが交絡されているか |
| 主な結果 | ブロック化された因果順序 + `NA` ブロックを含む隣接行列 | 祖先集合 + 交絡ペアが `NA` の隣接行列 |
| 答えられる問い | 「順序のどの部分なら信頼できるか」 | 「隠れた共通原因を共有しているのは正確にはどのペアか」 |

どちらもHSIC独立性検定に依存しており、そのコストは $`O(n^2)`$
で増加する。$`n`$ が
大きい場合はサブサンプリングを検討すること。関係が非線形なら
[CAM-UV](https://morimotoosamu.github.io/lingamr/articles/nonlinear-ja.md)
を使う。

## 関連記事

- [手法選択ガイド](https://morimotoosamu.github.io/lingamr/articles/method-selection-ja.md)
  — データに合う手法の選び方
- [非線形の手法（RESITとCAM-UV）](https://morimotoosamu.github.io/lingamr/articles/nonlinear-ja.md)
  — CAM-UVは非線形モデルの 潜在交絡変数を扱う
- [ブートストラップと診断](https://morimotoosamu.github.io/lingamr/articles/bootstrap-diagnostics-ja.md)
