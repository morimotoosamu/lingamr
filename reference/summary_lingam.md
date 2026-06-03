# Direct LiNGAM モデルの適合度を一括要約

推定済みの Direct LiNGAM モデルについて、LiNGAM
が依拠する2つの主要な前提
（残差の相互独立性・残差の非ガウス性）の成立度合いを一度に検証し、まとめて
表示する。内部で
[`get_error_independence_p_values()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values.md)
と
[`test_residual_normality()`](https://morimotoosamu.github.io/lingamr/reference/test_residual_normality.md)
を呼び出す。

## Usage

``` r
summary_lingam(
  X,
  lingam_result,
  independence_method = "spearman",
  normality_method = "shapiro",
  alpha = 0.05
)
```

## Arguments

- X:

  元データ (matrix or data.frame)。`lingam_result` の推定に用いたもの。

- lingam_result:

  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
  の返り値（`LingamResult` オブジェクト）

- independence_method:

  残差の独立性検定で用いる相関係数の種類 ("spearman", "pearson",
  "kendall")。[`get_error_independence_p_values()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values.md)
  に渡す。

- normality_method:

  残差の正規性検定の手法 ("shapiro", "ks", "ad", "lillie",
  "jb")。[`test_residual_normality()`](https://morimotoosamu.github.io/lingamr/reference/test_residual_normality.md)
  に渡す。

- alpha:

  有意水準 (default: 0.05)

## Value

`lingam_summary` クラスのリスト。以下の要素を含む：

- `n_variables`, `n_samples`: 変数の数・観測数

- `causal_order`: 因果順序（変数名ラベル）

- `n_edges`: 隣接行列の非ゼロ要素数（推定されたエッジ数）

- `independence_p_values`: 残差間の独立性検定の p 値行列

- `n_dependent_pairs`, `n_pairs`: p \< alpha のペア数 / 全ペア数

- `min_independence_p`: 独立性検定 p 値の最小値

- `normality`: 正規性検定の結果（`lingam_normality_test` オブジェクト）

- `n_non_gaussian`: 非ガウスと判定された変数の数

- `alpha`, `independence_method`, `normality_method`: 用いた設定

## Details

BIC/AIC のようなガウス尤度ベースの指標は、LiNGAM
の「誤差は非ガウス」という
前提と理論的に整合しないため含めていない。代わりに前提そのものの検証結果を
要約する。

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

model <- lingam_direct(LiNGAM_sample_1000$data, reg_method = "ols")

summary_lingam(LiNGAM_sample_1000$data, model)
#> === Direct LiNGAM Model Summary ===
#> Variables:    6
#> Observations: 1000
#> Edges:        15
#> Causal order: x3 -> x2 -> x0 -> x4 -> x5 -> x1
#> 
#> --- Assumption 1: Independence of residuals ---
#> Method:           spearman
#> Dependent pairs:  0 / 15  (p < 0.050)
#> Min p-value:      0.9187
#> => Residuals appear mutually independent (assumption supported).
#> 
#> --- Assumption 2: Non-Gaussianity of residuals ---
#> Method:           shapiro
#> Non-Gaussian:     6 / 6  (p <= 0.050)
#> => All residuals are non-Gaussian (assumption supported).
```
