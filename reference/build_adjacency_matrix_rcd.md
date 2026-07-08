# Build the RCD adjacency matrix from parents and confounder pairs

Faithful port of `rcd.py`'s adjacency-matrix construction (original
368-408 lines). `B[i, j]` is the coefficient of `j -> i`, matching the
lingamr convention; no transpose is needed.

## Usage

``` r
build_adjacency_matrix_rcd(X, P, C)
```

## Arguments

- X:

  (uncentered) data matrix

- P:

  parent list from
  [`extract_parents()`](https://morimotoosamu.github.io/lingamr/reference/extract_parents.md)

- C:

  confounder-pair list from
  [`extract_vars_sharing_confounders()`](https://morimotoosamu.github.io/lingamr/reference/extract_vars_sharing_confounders.md)

## Value

adjacency matrix B (n_features x n_features), with `NA` entries for
confounder pairs
