# LASSO regression (lambda selection by information criterion or CV)

LASSO regression (lambda selection by information criterion or CV)

## Usage

``` r
fit_lasso(y, Xp, lambda = "BIC")
```

## Arguments

- y:

  response variable

- Xp:

  predictor matrix

- lambda:

  lambda selection method "lambda.min" : minimum CV prediction error
  "lambda.1se" : CV 1SE rule "AIC" : minimum AIC "BIC" : minimum BIC,
  default

## Value

coefficient vector
