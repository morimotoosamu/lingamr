# Estimate the adjacency matrix from a causal order

Estimate the adjacency matrix from a causal order

## Usage

``` r
estimate_adjacency_matrix(
  X,
  causal_order,
  prior_knowledge = NULL,
  method = "adaptive_lasso",
  lambda = "BIC",
  init_method = "ols"
)
```

## Arguments

- X:

  original data

- causal_order:

  causal order (vector of 1-based indices)

- prior_knowledge:

  prior-knowledge matrix (NULL allowed)

- method:

  regression method "ols" : ordinary least squares (default) "lasso" :
  LASSO regression (glmnet) "adaptive_lasso": Adaptive LASSO (two-stage)
  "ridge" : Ridge regression (glmnet)

- lambda:

  LASSO penalty (NULL = automatic selection by cross-validation)
  "lambda.min" : minimum prediction error "lambda.1se" : 1SE rule
  (sparser) "AIC" : minimum AIC (no CV, fast) "BIC" : minimum BIC (no
  CV, fast, sparsest), default "oracle" : Adaptive LASSO only. Not
  usable with Ridge.

- init_method:

  estimation method for the initial weights of Adaptive LASSO "ols" :
  ordinary least squares (default) "ridge" : Ridge regression

## Value

adjacency matrix B (n_features x n_features)
