# DAG中の全パスを探索（from_index → to_index）

DAG中の全パスを探索（from_index → to_index）

## Usage

``` r
find_all_paths(adjacency_matrix, from_index, to_index, min_causal_effect = 0)
```

## Arguments

- adjacency_matrix:

  隣接行列

- from_index:

  始点インデックス (1-based)

- to_index:

  終点インデックス (1-based)

- min_causal_effect:

  最小因果効果の閾値

## Value

list(paths, effects)
