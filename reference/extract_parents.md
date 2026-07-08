# Extract parents from the ancestor sets

Faithful port of `rcd.py`'s `extract_parents()` (original 318-343
lines).

## Usage

``` r
extract_parents(X, M, cor_alpha)
```

## Arguments

- X:

  (uncentered) data matrix

- M:

  ancestor-set list from
  [`extract_ancestors()`](https://morimotoosamu.github.io/lingamr/reference/extract_ancestors.md)

- cor_alpha:

  correlation-test significance level

## Value

list of length `ncol(X)`; each element is a sorted integer vector of
parent indices (subset of the corresponding ancestor set)
