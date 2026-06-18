# VAR-LiNGAM for time series causal discovery

Fits a vector autoregressive (VAR) model to time series data and applies
Direct LiNGAM to the residuals to recover the instantaneous (lag-0)
causal structure. The lagged causal matrices are then derived from the
VAR coefficients and the instantaneous structure.

## Usage

``` r
lingam_var(
  X,
  lags = 1L,
  criterion = "bic",
  measure = "pwling",
  reg_method = "adaptive_lasso",
  lambda = "BIC",
  init_method = "ols",
  prune = TRUE
)
```

## Arguments

- X:

  numeric matrix or data frame (n_samples x n_features). Rows are
  ordered in time (earliest first).

- lags:

  maximum lag order. When `criterion` is not NULL, the best lag in
  `1:lags` is selected by the information criterion; otherwise `lags` is
  used directly.

- criterion:

  lag-selection criterion ("bic", "aic", "hqic", or "fpe"), or NULL to
  use `lags` directly without selection.

- measure:

  independence measure passed to
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
  ("pwling" or "kernel").

- reg_method:

  regression method for the instantaneous adjacency matrix:
  "adaptive_lasso" (default), "lasso", "ols", or "ridge" (see
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)).

- lambda:

  penalty (lambda) selection for the instantaneous matrix: "BIC"
  (default), "AIC", "lambda.min", "lambda.1se", or "oracle" (see
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)).

- init_method:

  initial-weight method for adaptive LASSO (see
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)).

- prune:

  logical; if `TRUE` (default, matching the Python reference), all
  adjacency matrices (instantaneous B0 and the lagged B_k) are refined
  together by adaptive LASSO so weak edges are shrunk toward zero.
  Requires the glmnet package. Set `FALSE` to keep the raw
  `B_k = (I - B0) M_k` matrices (no glmnet needed when
  `reg_method = "ols"`).

## Value

A `VARLiNGAMResult` object (list) containing:

- `adjacency_matrices`: array (1 + lags, n_features, n_features). The
  first slice `[1, , ]` is the instantaneous matrix B0; slice
  `[k + 1, , ]` is the lagged matrix B_k for lag k (k = 1..lags).
  Convention: `B[i, j]` is the effect from variable j to variable i.

- `causal_order`: estimated causal order of the instantaneous structure
  (1-based indices).

- `residuals`: VAR residuals (n_samples - lags, n_features).

- `lags`: the lag order actually used.

## Details

The model is `X_t = B0 X_t + sum_{k=1}^{p} B_k X_{t-k} + e_t`, where B0
is the instantaneous effect matrix (strictly acyclic) and e_t are
mutually independent non-Gaussian errors. VAR coefficients `M_k` are
estimated by ordinary least squares (no intercept); residuals
`e_t = X_t - sum M_k X_{t-k}` are passed to
[`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
to obtain B0, and the lagged matrices follow `B_k = (I - B0) M_k`.

## References

Hyvärinen, A., Zhang, K., Shimizu, S., & Hoyer, P. O. (2010). Estimation
of a structural vector autoregression model using non-Gaussianity.
*Journal of Machine Learning Research*, 11, 1709-1731. Ported from the
Python implementation cdt15/lingam (<https://github.com/cdt15/lingam>).
See also the VARLiNGAM R code of Moneta et al.
(<https://sites.google.com/site/dorisentner/publications/VARLiNGAM>).

## Examples

``` r
sample <- generate_varlingam_sample(n = 500, seed = 42)

# OLS instantaneous structure without pruning (no extra packages required)
model <- lingam_var(sample$data, lags = 1, reg_method = "ols", prune = FALSE)
round(model$adjacency_matrices[1, , ], 2)  # instantaneous B0
#>      x0   x1 x2
#> x0 0.00  0.0  0
#> x1 0.63  0.0  0
#> x2 0.02 -0.5  0
```
