# Compute all pairwise total causal effects at once

Equivalent to calling
[`calculate_total_effect()`](https://morimotoosamu.github.io/lingamr/reference/calculate_total_effect.md)
for every (from, to) pair, but computed as the truncated Neumann series
`B + B^2 + ... + B^(p-1)`. In a DAG every walk is a simple path (a
repeated node would imply a cycle), so the series contains exactly the
same path-product terms the DFS enumerates, and entries with no
connecting path are exactly zero. Only valid for acyclic `B` – with
cycles the series would count non-simple walks that
[`find_all_paths()`](https://morimotoosamu.github.io/lingamr/reference/find_all_paths.md)
excludes.

## Usage

``` r
calculate_total_effects_all(adjacency_matrix)
```

## Arguments

- adjacency_matrix:

  Adjacency matrix (p x p). `B[i,j]` is the coefficient of j -\> i; NAs
  are treated as absent edges.

## Value

p x p matrix `TE` with `TE[to, from]` = total effect of from on to
