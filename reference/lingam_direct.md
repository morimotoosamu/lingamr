# Direct LiNGAM

Direct LiNGAM

## Usage

``` r
lingam_direct(
  X,
  prior_knowledge = NULL,
  apply_prior_knowledge_softly = FALSE,
  measure = "pwling",
  reg_method = "adaptive_lasso",
  lambda = "BIC",
  init_method = "ols"
)
```

## Arguments

- X:

  数値行列 (n_samples x n_features), data frame or matrix

- prior_knowledge:

  事前知識行列 (n_features x n_features) または NULL。 0: x_i から x_j
  への有向パスなし 1: x_i から x_j への有向パスあり -1: 不明

- apply_prior_knowledge_softly:

  事前知識をソフトに適用するか (logical)

- measure:

  独立性の評価尺度 ("pwling" または "kernel")

- reg_method:

  隣接行列推定の回帰手法。 "ols": 最小二乗法、 "lasso": LASSO回帰、
  "adaptive_lasso": 適応的LASSO回帰（デフォルト）。

- lambda:

  LASSO のペナルティ（ラムダ）選択。 "lambda.min" : CV予測誤差最小,
  予測精度優先。 "lambda.1se" : CV 1SEルール、ロバストで過学習しにくい。
  "AIC": AIC最小。高速。 "BIC":
  BIC最小。高速、最もスパース。デフォルト。 "oracle"
  ：適応的LASSO回帰のみ。オラクル性を担保したλを選択。高速。

- init_method:

  適応的LASSO回帰の初期重みの推定手法。 "ols":
  最小二乗法（デフォルト）、 "ridge": Ridege回帰。
  多重共線性が疑われる場合はRidege回帰がおすすめ。

## Value

list(adjacency_matrix, causal_order)

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

# OLS (no additional packages required)
result <- lingam_direct(LiNGAM_sample_1000$data, reg_method = "ols")
round(result$adjacency_matrix, 3)
#>       x0 x1     x2     x3     x4    x5
#> x0 0.000  0 -0.040  3.274  0.000 0.000
#> x1 3.237  0  1.965  0.014 -0.034 0.006
#> x2 0.000  0  0.000  5.993  0.000 0.000
#> x3 0.000  0  0.000  0.000  0.000 0.000
#> x4 7.992  0 -1.062  0.394  0.000 0.000
#> x5 3.873  0  0.069 -0.315  0.018 0.000

# \donttest{
# LASSO (requires glmnet)
result_lasso <- lingam_direct(LiNGAM_sample_1000$data)
round(result_lasso$adjacency_matrix, 3)
#>       x0 x1     x2    x3 x4 x5
#> x0 0.000  0  0.000 3.033  0  0
#> x1 2.988  0  2.002 0.000  0  0
#> x2 0.000  0  0.000 5.993  0  0
#> x3 0.000  0  0.000 0.000  0  0
#> x4 8.017  0 -1.009 0.000  0  0
#> x5 4.015  0  0.000 0.000  0  0
# }
```
