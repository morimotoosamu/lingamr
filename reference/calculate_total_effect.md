# 隣接行列から2変数間の総合因果効果を計算する

[`find_all_paths()`](https://morimotoosamu.github.io/lingamr/reference/find_all_paths.md)
で列挙した全経路効果の総和を返す。 パスが存在しない場合は 0 を返す。

## Usage

``` r
calculate_total_effect(adjacency_matrix, from_index, to_index)
```

## Arguments

- adjacency_matrix:

  隣接行列 (n x n)。`B[i,j]` は j → i の係数。

- from_index:

  原因変数のインデックス (1-based)

- to_index:

  結果変数のインデックス (1-based)

## Value

総合因果効果（スカラー）
