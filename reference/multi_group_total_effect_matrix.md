# Compute a total-effect matrix from an adjacency matrix via path products

Unlike
[`estimate_all_total_effects()`](https://morimotoosamu.github.io/lingamr/reference/estimate_all_total_effects.md)
(regression-based), this sums path products over the DAG defined by `B`,
matching the upstream Python MultiGroup bootstrap's
[`calculate_total_effect()`](https://morimotoosamu.github.io/lingamr/reference/calculate_total_effect.md).

## Usage

``` r
multi_group_total_effect_matrix(B, causal_order)
```

## Arguments

- B:

  adjacency matrix (n_features x n_features), `B[i,j]` = j -\> i

- causal_order:

  causal order (1-based indices)

## Value

total-effect matrix (n_features x n_features)
