# Fit a VAR(p) model by OLS (no intercept)

Fit a VAR(p) model by OLS (no intercept)

## Usage

``` r
fit_var_ols(X, lags)
```

## Arguments

- X:

  numeric matrix (n_samples x n_features), rows ordered in time

- lags:

  lag order (positive integer)

## Value

list with `coefs` (array (lags, n_features, n_features); `coefs[k, , ]`
is M_k such that `X_t = sum_k M_k X_{t-k} + e_t`) and `residuals`
(n_samples - lags, n_features)
