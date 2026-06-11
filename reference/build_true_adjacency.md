# エッジ指定から真の隣接行列を構築する

エッジ指定から真の隣接行列を構築する

## Usage

``` r
build_true_adjacency(var_names, from, to, coef)
```

## Arguments

- var_names:

  変数名ベクトル

- from:

  エッジの原因変数名ベクトル

- to:

  エッジの結果変数名ベクトル（from と同じ長さ）

- coef:

  エッジ係数ベクトル（from と同じ長さ）

## Value

隣接行列 (p x p)。`m[to, from] = coef`（行 = to, 列 = from）
