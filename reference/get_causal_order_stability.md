# ブートストラップによる因果順序の安定性を評価

各ブートストラップ標本で推定された因果順序 (causal_order)
を集計し、順序が
どれだけ安定しているかを数値化する。各変数の順位分布、変数ペアの先行確率、
および全体の安定性スコアを返す。

## Usage

``` r
get_causal_order_stability(result, labels = NULL)
```

## Arguments

- result:

  BootstrapResult オブジェクト（現行バージョンで実行したもの）

- labels:

  変数名ベクトル (NULL の場合は x0, x1, ... を自動生成)

## Value

`causal_order_stability` クラスのリスト。以下を含む：

- `rank_summary`: 各変数の順位の要約 (variable, mean_rank, sd_rank,
  median_rank, mode_rank)。mean_rank 昇順（上流から）にソート済み。
  順位は 1 が最も上流。

- `precedence_matrix`: 先行確率行列。`P[i, j]` は変数 i が変数 j
  より上流 （先）に位置したブートストラップ標本の割合。

- `stability_score`: 全体の安定性スコア。0（順序がランダム）〜
  1（全標本で順序が一致）。各変数ペアの先行確率が 0/1 に近いほど高い。

- `n_sampling`: ブートストラップ標本数。

## Examples

``` r
dat <- generate_lingam_sample_6()
bs <- lingam_direct_bootstrap(dat$data, n_sampling = 30L, seed = 42)
#> Bootstrap: 30 iterations, method=adaptive_lasso (sequential)
#>   iteration 1 / 30
#>   iteration 10 / 30
#>   iteration 20 / 30
#>   iteration 30 / 30
#> Completed in 0.9 seconds.
get_causal_order_stability(bs, labels = names(dat$data))
#> === Causal Order Stability ===
#> Bootstrap samples:       30
#> Overall stability score: 0.680  (0 = random, 1 = fully stable)
#> 
#> Rank summary (sorted by mean rank; 1 = most upstream):
#>  variable mean_rank sd_rank median_rank mode_rank
#>        x3      1.17    0.91         1.0         1
#>        x0      2.57    0.57         3.0         3
#>        x2      2.93    0.98         2.5         2
#>        x5      4.33    1.32         4.0         3
#>        x4      4.87    0.86         5.0         5
#>        x1      5.13    1.11         5.0         6
#> 
#> Precedence probability P[i, j] = P(variable i precedes j):
#>      x0   x1   x2   x3   x4   x5
#> x0 0.00 0.97 0.47 0.03 0.97 1.00
#> x1 0.03 0.00 0.03 0.03 0.40 0.37
#> x2 0.53 0.97 0.00 0.03 0.97 0.57
#> x3 0.97 0.97 0.97 0.00 0.97 0.97
#> x4 0.03 0.60 0.03 0.03 0.00 0.43
#> x5 0.00 0.63 0.43 0.03 0.57 0.00
```
