# Total causal effect from a joined adjacency matrix (graph-based)

Computes the total effect by summing path products over the
time-expanded graph, reusing
[`calculate_total_effect()`](https://morimotoosamu.github.io/lingamr/reference/calculate_total_effect.md).
Port of the Python reference `estimate_total_effect2`; used internally
by the VAR-LiNGAM bootstrap.

## Usage

``` r
var_total_effect_graph(am_joined, from_index, to_index)
```

## Arguments

- am_joined:

  joined adjacency matrix (n_features x n_features\*(1 + lags)),
  `B[i, j]` is the coefficient of j -\> i

- from_index:

  source column in the joined index space (1-based)

- to_index:

  destination column in the joined index space (1-based)

## Value

the total effect (scalar)
