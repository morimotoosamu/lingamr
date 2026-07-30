# Prune VARMA-LiNGAM adjacency matrices by adaptive LASSO

Re-estimates the instantaneous matrix B0, every lagged matrix psi_tau,
and every MA matrix omega_w jointly, shrinking weak edges to zero. Port
of the Python reference `_pruning`. For each target variable, the
predictors are its contemporaneous ancestors (those preceding it in
`causal_order`) plus all variables at lags 1..p and the LiNGAM residuals
e_t at lags 1..q. Unlike the Python reference (which wraps rows around
with `np.roll`), the regression window is restricted to rows where every
regressor is observed.

## Usage

``` r
prune_varma_lingam(
  X,
  ee_full,
  causal_order,
  order,
  lambda = "BIC",
  init_method = "ols"
)
```

## Arguments

- X:

  numeric matrix (n_samples x n_features), rows ordered in time

- ee_full:

  LiNGAM residuals `e_t = (I - B0) n_t`, full length (n_samples rows;
  the first max(p, q) rows are zero)

- causal_order:

  instantaneous causal order (1-based indices)

- order:

  VARMA order c(p, q)

- lambda:

  lambda selection passed to
  [`fit_adaptive_lasso()`](https://morimotoosamu.github.io/lingamr/reference/fit_adaptive_lasso.md)

- init_method:

  initial-weight method for adaptive LASSO

## Value

list with `psis` (array (1 + p, m, m); slice 1 is B0) and `omegas`
(array (q, m, m))
