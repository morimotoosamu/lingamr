# 因果順序から隣接行列を推定

因果順序から隣接行列を推定

## Usage

``` r
estimate_adjacency_matrix(
  X,
  causal_order,
  prior_knowledge = NULL,
  method = "adaptive_lasso",
  lambda = "BIC",
  init_method = "ols"
)
```

## Arguments

- X:

  元データ

- causal_order:

  因果順序 (1-based index のベクトル)

- prior_knowledge:

  事前知識行列 (NULL可)

- method:

  回帰手法 "ols" : 通常の最小二乗法（デフォルト） "lasso" :
  LASSO回帰（glmnet） "adaptive_lasso": Adaptive LASSO（2段階）

- lambda:

  LASSO のペナルティ (NULL = 交差検証で自動選択) "lambda.min" :
  予測誤差最小 "lambda.1se" : 1SE ルール（よりスパース） "AIC" :
  AIC最小（CVなし、高速） "BIC" :
  BIC最小（CVなし、高速、最もスパース）デフォルト

- init_method:

  Adaptive LASSOの初期重みの推定手法 "ols" :最小二乗法（デフォルト）
  "ridge" :Ridge回帰

## Value

隣接行列 B (n_features x n_features)
