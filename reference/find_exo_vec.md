# Find the most sink-like candidate variable

Faithful port of `_search_causal_order.find_exo_vec()` (original 227-276
lines). For each candidate `j` in `Uc`, regresses `j` on the other
candidates (`setdiff(Uc, j)`, **not** `U`) and evaluates how independent
the residual is from those explanatory variables.

## Usage

``` r
find_exo_vec(X, Uc, U, independence, Cov)
```

## Arguments

- X:

  data matrix

- Uc:

  candidate variable indices

- U:

  all currently undetermined variable indices

- independence:

  "hsic" or "fcorr"

- Cov:

  precomputed `stats::cov(X)`, since `X` is invariant across all calls
  within a single
  [`parce_search_causal_order()`](https://morimotoosamu.github.io/lingamr/reference/parce_search_causal_order.md)
  search

## Value

list(m = selected variable index, eval = its evaluation value)
