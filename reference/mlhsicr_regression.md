# HSIC-sum-minimizing ("MLHSICR") regression

Regresses `Y[, xi]` on `Y[, xj_list]` (no intercept) by minimizing the
sum of the empirical HSIC statistics between the residual and each
explanatory variable, instead of ordinary least squares. Used as a
fallback when the OLS residual is not independent of the explanatory
variables.

## Usage

``` r
mlhsicr_regression(Y, xi, xj_list)
```

## Arguments

- Y:

  data matrix (residualized already, if applicable)

- xi:

  target column index

- xj_list:

  explanatory column indices (length \>= 2)

## Value

list(resid = residual vector, coef = coefficient vector)

## Details

The kernel width used to build the residual's Gram matrix is itself a
linear combination of the explanatory variables' kernel widths (faithful
to the upstream implementation; see original 207-208 lines). This is an
unusual design choice but is deliberately preserved as-is.
