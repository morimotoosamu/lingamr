# Adaptive LASSO

Adaptive LASSO

## Usage

``` r
fit_adaptive_lasso(
  y,
  Xp,
  lambda = "BIC",
  gamma_weight = 1,
  init_method = "ols"
)
```

## Arguments

- y:

  response variable

- Xp:

  predictor matrix

- lambda:

  lambda selection method ("lambda.min", "lambda.1se", "AIC", "BIC",
  "oracle")

- gamma_weight:

  exponent of the weights

- init_method:

  estimation method for the initial weights ("ols" or "ridge")

## Value

coefficient vector
