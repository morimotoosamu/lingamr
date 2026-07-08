# Causal order search via pwling, jointly across multiple groups

Sample-size-weighted sum of the single-group pwling objective
([`search_causal_order_pwling()`](https://morimotoosamu.github.io/lingamr/reference/search_causal_order_pwling.md))
across all groups in `X_list`, following the joint estimation objective
of Shimizu (2012). Reuses the same standardize-once / correlation-matrix
/ antisymmetry optimizations as
[`search_causal_order_pwling()`](https://morimotoosamu.github.io/lingamr/reference/search_causal_order_pwling.md),
applied independently within each group before the weighted sum.

## Usage

``` r
search_causal_order_pwling_multi(X_list, U, Uc, Vj)
```

## Arguments

- X_list:

  List of per-group data matrices (residualized so far), one per group.
  All must have the same number of columns.

- U:

  Indices of all currently undetermined variables (shared across groups)

- Uc:

  Indices of candidate variables (shared across groups)

- Vj:

  Variable set based on prior knowledge (shared across groups)

## Value

Index of the selected variable
