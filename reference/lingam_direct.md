# Direct LiNGAM

Direct LiNGAM

## Usage

``` r
lingam_direct(
  X,
  prior_knowledge = NULL,
  apply_prior_knowledge_softly = FALSE,
  measure = "pwling",
  reg_method = "adaptive_lasso",
  lambda = "BIC",
  init_method = "ols"
)
```

## Arguments

- X:

  Numeric matrix (n_samples x n_features), data frame or matrix

- prior_knowledge:

  Prior knowledge matrix (n_features x n_features) or NULL. 0: no
  directed path from x_i to x_j 1: directed path from x_i to x_j -1:
  unknown

- apply_prior_knowledge_softly:

  Whether to apply prior knowledge softly (logical)

- measure:

  Independence evaluation measure ("pwling" or "kernel")

- reg_method:

  Regression method for adjacency matrix estimation. "ols": ordinary
  least squares, "lasso": LASSO regression, "adaptive_lasso": adaptive
  LASSO regression (default), "ridge": Ridge regression (robust to
  multicollinearity; does not perform sparse estimation).

- lambda:

  LASSO penalty (lambda) selection. "lambda.min" : minimum CV prediction
  error, prioritizes prediction accuracy. "lambda.1se" : CV 1SE rule,
  robust and less prone to overfitting. "AIC": minimum AIC. Fast. "BIC":
  minimum BIC. Fast, sparsest. Default. "oracle" : adaptive LASSO
  regression only. Selects a lambda that guarantees the oracle property.
  Fast.

- init_method:

  Method for estimating the initial weights of adaptive LASSO
  regression. "ols": ordinary least squares (default), "ridge": Ridge
  regression. Ridge regression is recommended when multicollinearity is
  suspected.

## Value

A `LingamResult` object (list) containing the following elements:

- `adjacency_matrix`: adjacency matrix B (n_features x n_features).
  **Convention: `B[i, j]` is the causal coefficient from variable j to
  variable i (j -\> i).** Zero elements indicate no causal relationship.

- `causal_order`: estimated causal order (integer vector of 1-based
  indices). Earlier elements are more upstream (closer to exogenous
  variables).

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

# OLS (no additional packages required)
result <- lingam_direct(LiNGAM_sample_1000$data, reg_method = "ols")
round(result$adjacency_matrix, 3)
#>       x0 x1     x2     x3     x4    x5
#> x0 0.000  0 -0.040  3.274  0.000 0.000
#> x1 3.237  0  1.965  0.014 -0.034 0.006
#> x2 0.000  0  0.000  5.993  0.000 0.000
#> x3 0.000  0  0.000  0.000  0.000 0.000
#> x4 7.992  0 -1.062  0.394  0.000 0.000
#> x5 3.873  0  0.069 -0.315  0.018 0.000

# \donttest{
# LASSO (requires glmnet)
result_lasso <- lingam_direct(LiNGAM_sample_1000$data)
round(result_lasso$adjacency_matrix, 3)
#>       x0 x1     x2    x3 x4 x5
#> x0 0.000  0  0.000 3.033  0  0
#> x1 2.988  0  2.002 0.000  0  0
#> x2 0.000  0  0.000 5.993  0  0
#> x3 0.000  0  0.000 0.000  0  0
#> x4 8.017  0 -1.009 0.000  0  0
#> x5 4.015  0  0.000 0.000  0  0
# }
```
