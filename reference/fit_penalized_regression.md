# Penalized regression via glmnet (IC or CV lambda selection)

Internal helper shared by
[`fit_lasso()`](https://morimotoosamu.github.io/lingamr/reference/fit_lasso.md)
and
[`fit_ridge_reg()`](https://morimotoosamu.github.io/lingamr/reference/fit_ridge_reg.md).
Both functions differ only in `alpha` and `lambda_seq`; this function
encapsulates the duplicated IC / CV branches.

## Usage

``` r
fit_penalized_regression(y, Xp_mat, alpha, lambda, lambda_seq)
```

## Arguments

- y:

  response variable (numeric vector)

- Xp_mat:

  predictor matrix (already coerced to matrix)

- alpha:

  glmnet mixing parameter: 1 = LASSO, 0 = Ridge

- lambda:

  lambda selection method ("AIC", "BIC", "lambda.min", "lambda.1se")

- lambda_seq:

  numeric vector of lambda values passed to glmnet

## Value

coefficient vector (excluding intercept)
