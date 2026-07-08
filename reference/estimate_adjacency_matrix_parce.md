# Estimate the adjacency matrix from a ParceLiNGAM causal order

`causal_order` is a list whose first element may be an unresolved block
(length \> 1); all remaining elements are length-1. Block members are
never regression targets (their parents cannot be identified), but they
are valid predictors for downstream variables. Pairs within the block
are set to `NA`.

## Usage

``` r
estimate_adjacency_matrix_parce(
  X,
  causal_order,
  prior_knowledge,
  method,
  lambda,
  init_method
)
```

## Arguments

- X:

  original (uncentered) data

- causal_order:

  list as produced by
  [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)

- prior_knowledge:

  prior-knowledge matrix (NULL allowed)

- method:

  regression method

- lambda:

  lambda selection

- init_method:

  adaptive LASSO initial-weight method

## Value

adjacency matrix B (n_features x n_features)
