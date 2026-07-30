# Residual of one variable regressed on a set of explanatory variables

Faithful port of `CAMUV._get_residual()`. With no explanatory variables
the column is returned as is.

## Usage

``` r
camuv_get_residual(X, explained_i, explanatory_ids, reg_fn)
```

## Arguments

- X:

  data matrix

- explained_i:

  column index of the explained variable

- explanatory_ids:

  column indices of the explanatory variables (possibly empty)

- reg_fn:

  regressor function from
  [`resit_make_regressor()`](https://morimotoosamu.github.io/lingamr/reference/resit_make_regressor.md)

## Value

residual vector (length n)
