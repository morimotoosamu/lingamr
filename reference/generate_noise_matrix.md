# 変数ごとに独立シードでノイズ行列を生成する

列 k は `set.seed(seed + k - 1)` の直後に `noise_fn(n)` で生成される。
変数ごとにシードを固定することで、同じ seed
なら常に同じノイズ列が得られる。

## Usage

``` r
generate_noise_matrix(n, n_vars, seed, noise_fn)
```

## Arguments

- n:

  サンプルサイズ

- n_vars:

  変数（列）の数

- seed:

  基準シード。列 k には seed + k - 1 を用いる

- noise_fn:

  `function(n)` 形式のノイズ生成関数

## Value

ノイズ行列 (n x n_vars)
