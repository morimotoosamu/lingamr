# Adaptive LASSO

Adaptive LASSO

## Usage

``` r
fit_adaptive_lasso(
  y,
  Xp,
  lambda = "BIC",
  gamma_weight = 1,
  init_method = "ols"
)
```

## Arguments

- y:

  目的変数

- Xp:

  説明変数行列

- lambda:

  ラムダ選択方法 ("lambda.min", "lambda.1se", "AIC", "BIC", "oracle")

- gamma_weight:

  重みの指数

- init_method:

  初期重みの推定手法 ("ols" または "ridge")

## Value

係数ベクトル
