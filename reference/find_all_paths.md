# DAG 中の全パスを深さ優先探索で列挙する

`B[i, j]` が j → i を表す隣接行列を受け取り、`from_index` から
`to_index` に至る全パスとそれぞれの経路効果（係数の積）を返す。

## Usage

``` r
find_all_paths(adjacency_matrix, from_index, to_index, min_causal_effect = 0)
```

## Arguments

- adjacency_matrix:

  隣接行列 (n x n)。`B[i,j]` は j → i の係数。

- from_index:

  始点インデックス (1-based)

- to_index:

  終点インデックス (1-based)

- min_causal_effect:

  このしきい値以下の係数は存在しないエッジとみなす

## Value

list(paths, effects)
