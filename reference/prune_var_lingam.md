# Prune VAR-LiNGAM adjacency matrices by adaptive LASSO

Re-estimates the instantaneous matrix B0 and every lagged matrix B_k
jointly, shrinking weak edges to zero. Port of the Python reference
`_pruning`. For each target variable, the predictors are its
contemporaneous ancestors (those preceding it in `causal_order`) plus
all variables at lags 1..lags; the coefficients are fitted by adaptive
LASSO and written back into B.

## Usage

``` r
prune_var_lingam(X, causal_order, lags, lambda = "BIC", init_method = "ols")
```

## Arguments

- X:

  numeric matrix (n_samples x n_features), rows ordered in time

- causal_order:

  instantaneous causal order (1-based indices)

- lags:

  lag order

- lambda:

  lambda selection passed to
  [`fit_adaptive_lasso()`](https://morimotoosamu.github.io/lingamr/reference/fit_adaptive_lasso.md)

- init_method:

  initial-weight method for adaptive LASSO

## Value

array (lags + 1, n_features, n_features); slice 1 is B0, slice k+1 is
B_k
