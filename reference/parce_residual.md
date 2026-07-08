# Residual of `X[, j]` regressed on `X[, xi_index]` via the pseudo-inverse of the covariance matrix

Residual of `X[, j]` regressed on `X[, xi_index]` via the pseudo-inverse
of the covariance matrix

## Usage

``` r
parce_residual(X, xi_index, j, Cov)
```

## Arguments

- X:

  data matrix

- xi_index:

  explanatory-variable indices (may be empty)

- j:

  target variable index

- Cov:

  precomputed `stats::cov(X)`, since `X` is invariant across all calls
  within a single
  [`parce_search_causal_order()`](https://morimotoosamu.github.io/lingamr/reference/parce_search_causal_order.md)
  search

## Value

residual vector
