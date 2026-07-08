# Adaptive LASSO with CV-selected lambda (n \<= p route)

Replicates upstream `_predict_adaptive_lasso` (StandardScaler + OLS
weights + `LassoLarsCV`), substituting `glmnet::cv.glmnet(alpha = 1)`
for `LassoLarsCV` (see
[`lingam_high_dim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_high_dim.md)
Details for the rationale).

## Usage

``` r
predict_adaptive_lasso_cv(X, predictors, target, gamma = 1)
```

## Arguments

- X:

  original-scale data matrix

- predictors:

  indices of predictor variables (1-based)

- target:

  index of the target variable (1-based)

- gamma:

  exponent of the adaptive weights (fixed at 1.0 upstream)

## Value

coefficient vector, same length and order as `predictors`
