# Scale factor used to make the (otherwise fixed, absolute) lambda search grids adapt to the response's natural scale.

`lasso_lambda_seq` / `ridge_lambda_seq` are fixed absolute grids.
Because glmnet's `standardize = TRUE` only standardizes the predictors
(not the response `y`), the penalty strength needed for meaningful
shrinkage scales with the magnitude of `y`. Without this scaling,
multiplying the whole input data by a constant changes which edges the
default `adaptive_lasso` + `lambda = "BIC"`/`"AIC"` pipeline selects,
even though the underlying relationships are identical up to that
constant.

## Usage

``` r
lambda_scale_factor(y)
```

## Arguments

- y:

  response variable (numeric vector)

## Value

a positive scale factor
