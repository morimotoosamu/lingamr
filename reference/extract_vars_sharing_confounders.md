# Detect variable pairs sharing an unobserved latent confounder

Faithful port of `rcd.py`'s `extract_vars_sharing_confounders()`
(original 352-366 lines). Only pairs with no parent-child relationship
in either direction are considered.

## Usage

``` r
extract_vars_sharing_confounders(X, P, cor_alpha)
```

## Arguments

- X:

  (uncentered) data matrix

- P:

  parent list from
  [`extract_parents()`](https://morimotoosamu.github.io/lingamr/reference/extract_parents.md)

- cor_alpha:

  correlation-test significance level

## Value

list of length `ncol(X)`; each element is a sorted integer vector of
indices sharing a latent confounder with that variable (symmetric)
