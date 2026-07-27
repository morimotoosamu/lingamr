# Dispatch a single regression to the backend selected by `method`

Central dispatcher shared by every place that fits "y on Xp with the
user-chosen regression method" (adjacency estimation, total effects,
Parce/RCD variants). Callers are expected to have validated `method`,
`lambda`, and `init_method` already; no validation happens here so that
error behaviour stays with the caller.

## Usage

``` r
fit_coef_by_method(y, Xp, method, lambda, init_method)
```

## Arguments

- y:

  response variable (numeric vector)

- Xp:

  predictor matrix

- method:

  one of "ols", "lasso", "adaptive_lasso", "ridge"

- lambda:

  lambda selection rule (ignored for OLS)

- init_method:

  initial estimator for adaptive LASSO (ignored otherwise)

## Value

coefficient vector (excluding intercept)
