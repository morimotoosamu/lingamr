# LASSO 回帰（情報量基準 or CV でラムダ選択）

LASSO 回帰（情報量基準 or CV でラムダ選択）

## Usage

``` r
fit_lasso(y, Xp, lambda = "BIC")
```

## Arguments

- y:

  目的変数

- Xp:

  説明変数行列

- lambda:

  ラムダ選択方法 "lambda.min" : CV予測誤差最小 "lambda.1se" : CV
  1SEルール "AIC" : AIC最小 "BIC" : BIC最小。デフォルト

## Value

係数ベクトル
