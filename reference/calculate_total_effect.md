# Compute the total causal effect between two variables from an adjacency matrix

Returns the sum of all path effects enumerated by
[`find_all_paths()`](https://morimotoosamu.github.io/lingamr/reference/find_all_paths.md).
Returns 0 if no path exists.

## Usage

``` r
calculate_total_effect(adjacency_matrix, from_index, to_index)
```

## Arguments

- adjacency_matrix:

  Adjacency matrix (n x n). `B[i,j]` is the coefficient of j -\> i.

- from_index:

  Index of the cause variable (1-based)

- to_index:

  Index of the effect variable (1-based)

## Value

Total causal effect (scalar)
