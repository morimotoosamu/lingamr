# Whether the OLS (or MLHSICR) residual of `xi` on `xj_list` is independent of every `xj`

Faithful port of `rcd.py`'s `is_independent()` helper (original 226-252
lines).

## Usage

``` r
is_independent_of_resid(
  Y,
  xi,
  xj_list,
  MLHSICR,
  independence,
  ind_alpha,
  ind_corr
)
```

## Arguments

- Y:

  residual matrix (from
  [`extract_ancestors()`](https://morimotoosamu.github.io/lingamr/reference/extract_ancestors.md))

- xi:

  candidate sink variable

- xj_list:

  explanatory-variable indices

- MLHSICR:

  whether to retry with
  [`mlhsicr_regression()`](https://morimotoosamu.github.io/lingamr/reference/mlhsicr_regression.md)
  on failure

- independence:

  "hsic" or "fcorr"

- ind_alpha:

  significance level (hsic only)

- ind_corr:

  rejection threshold (fcorr only)

## Value

TRUE if independent
