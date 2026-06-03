# 事前知識行列を作成

事前知識行列を作成

## Usage

``` r
make_prior_knowledge(
  n_variables,
  exogenous_variables = NULL,
  sink_variables = NULL,
  paths = NULL,
  no_paths = NULL,
  labels = NULL
)
```

## Arguments

- n_variables:

  変数の数

- exogenous_variables:

  外生変数 (1-based index または変数名, NULL可)
  指定した変数は他のどの変数からも影響を受けないとする

- sink_variables:

  シンク変数 (1-based index または変数名, NULL可)
  指定した変数は他のどの変数にも影響を与えないとする

- paths:

  有向パスが存在する変数ペア (NULL可) list(c(from, to), ...)
  の形式。インデックスまたは変数名で指定

- no_paths:

  有向パスが存在しない変数ペア (NULL可) list(c(from, to), ...)
  の形式。インデックスまたは変数名で指定

- labels:

  変数名ベクトル (NULL可) 変数名で指定する場合は必須。data.frame の
  colnames() 等を渡す

## Value

事前知識行列 (n_variables x n_variables) -1: 不明, 0: パスなし, 1:
パスあり

## Examples

``` r
# インデックスで指定
pk <- make_prior_knowledge(6, exogenous_variables = c(4))

# 変数名で指定
pk <- make_prior_knowledge(6,
  exogenous_variables = "x3",
  sink_variables = c("x1", "x4"),
  paths = list(c("x3", "x0"), c("x3", "x2")),
  no_paths = list(c("x5", "x2")),
  labels = c("x0", "x1", "x2", "x3", "x4", "x5")
)
```
