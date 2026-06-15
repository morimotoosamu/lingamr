# Enumerate all paths in a DAG via depth-first search

Takes an adjacency matrix where `B[i, j]` represents j -\> i, and
returns all paths from `from_index` to `to_index` together with each
path effect (the product of the coefficients).

## Usage

``` r
find_all_paths(adjacency_matrix, from_index, to_index, min_causal_effect = 0)
```

## Arguments

- adjacency_matrix:

  Adjacency matrix (n x n). `B[i,j]` is the coefficient of j -\> i.

- from_index:

  Start index (1-based)

- to_index:

  End index (1-based)

- min_causal_effect:

  Coefficients at or below this threshold are treated as nonexistent
  edges

## Value

list(paths, effects)
