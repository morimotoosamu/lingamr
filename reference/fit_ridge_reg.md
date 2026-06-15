# Ridge regression (lambda selection by information criterion or CV)

Ridge regression (lambda selection by information criterion or CV)

## Usage

``` r
fit_ridge_reg(y, Xp, lambda = "BIC")
```

## Arguments

- y:

  response variable

- Xp:

  predictor matrix

- lambda:

  lambda selection method "lambda.min" : minimum CV prediction error
  "lambda.1se" : CV 1SE rule "AIC" : minimum AIC "BIC" : minimum BIC,
  default "oracle" is not usable (Adaptive LASSO only).

## Value

coefficient vector
