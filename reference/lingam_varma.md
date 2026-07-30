# VARMA-LiNGAM for time series causal discovery

Fits a vector autoregressive moving-average (VARMA) model to time series
data and applies Direct LiNGAM to the residuals to recover the
instantaneous (lag-0) causal structure. The lagged causal matrices (psi)
and the moving-average causal matrices (omega) are then derived from the
VARMA coefficients and the instantaneous structure.

## Usage

``` r
lingam_varma(
  X,
  order = c(1L, 1L),
  criterion = "bic",
  measure = "pwling",
  reg_method = "adaptive_lasso",
  lambda = "BIC",
  init_method = "ols",
  prune = TRUE,
  ar_coefs = NULL,
  ma_coefs = NULL
)
```

## Arguments

- X:

  numeric matrix or data frame (n_samples x n_features). Rows are
  ordered in time (earliest first).

- order:

  VARMA order `c(p, q)` (AR and MA lags, non-negative integers, not both
  zero). When `criterion` is not NULL, the best order in `0:p x 0:q`
  (excluding `c(0, 0)`) is selected by the information criterion;
  otherwise `order` is used directly.

- criterion:

  order-selection criterion ("bic", "aic", or "hqic"), or NULL to use
  `order` directly without selection.

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
  adjacency matrices (psi and omega) are refined together by adaptive
  LASSO so weak edges are shrunk toward zero. Requires the glmnet
  package. Set `FALSE` to keep the raw transformed matrices (no glmnet
  needed when `reg_method = "ols"`).

- ar_coefs:

  optional known AR coefficient array (p, n_features, n_features) in
  reduced form. Must be supplied together with `ma_coefs`; when both are
  given the VARMA estimation step is skipped, `order` is derived from
  the arrays, and `criterion` is ignored.

- ma_coefs:

  optional known MA coefficient array (q, n_features, n_features) in
  reduced form. See `ar_coefs`.

## Value

A `VARMALiNGAMResult` object (list) containing:

- `adjacency_matrices`: list with two arrays. `psis` (1 + p, n_features,
  n_features): slice `[1, , ]` is the instantaneous matrix B0 (= psi_0)
  and slice `[tau + 1, , ]` is the lagged matrix psi_tau. `omegas` (q,
  n_features, n_features): slice `[w, , ]` is the MA causal matrix
  omega_w. Convention: entry `[i, j]` is the effect from variable j to
  variable i.

- `causal_order`: estimated causal order of the instantaneous structure
  (1-based indices).

- `residuals`: filtered VARMA residuals n_t (n_samples - max(p, q),
  n_features).

- `order`: the VARMA order `c(p, q)` actually used.

- `ar_coefs`, `ma_coefs`: reduced-form coefficient arrays Phi (p, m, m)
  and Theta (q, m, m) used to derive the causal matrices.

- `const`: estimated intercept vector of the reduced-form VARMA (zero
  when `ar_coefs`/`ma_coefs` are supplied).

## Details

The structural model is
`x_t = B0 x_t + sum_{tau=1}^{p} psi_tau x_{t-tau} + e_t + sum_{w=1}^{q} omega_w e_{t-w}`,
where B0 is the instantaneous effect matrix (strictly acyclic) and e_t
are mutually independent non-Gaussian errors. The reduced-form VARMA
`x_t = c + sum Phi_tau x_{t-tau} + n_t + sum Theta_w n_{t-w}` is
estimated first; Direct LiNGAM applied to the residuals n_t yields B0,
and the causal matrices follow `psi_tau = (I - B0) Phi_tau` and
`omega_w = (I - B0) Theta_w (I - B0)^{-1}`.

Differences from the Python reference implementation:

- VARMA coefficients are estimated by the two-stage Hannan-Rissanen
  procedure (a long VAR yields residual estimates, then X is regressed
  on its own lags and the lagged residuals, both stages with an
  intercept) instead of state-space maximum likelihood (statsmodels
  `VARMAX`). Estimates are consistent but not identical to the Python
  output; the `max_iter` parameter of the Python class has no
  counterpart here because the procedure is non-iterative.

- The MA causal matrices use the full similarity transform
  `omega_w = (I - B0) Theta_w (I - B0)^{-1}`. The Python implementation
  drops the trailing inverse factor due to a three-argument `np.dot`
  call (the third argument is treated as an output buffer), so its
  unpruned omegas equal `(I - B0) Theta_w`.

- Residual filtering initializes the first `max(p, q)` residuals with
  zeros instead of random draws, making the fit deterministic.

## References

Kawahara, Y., Shimizu, S., & Washio, T. (2011). Analyzing relationships
among ARMA processes based on non-Gaussianity of external influences.
*Neurocomputing*, 74(12-13), 2212-2221. Ported from the Python
implementation cdt15/lingam (<https://github.com/cdt15/lingam>).

## Examples

``` r
sample <- generate_varmalingam_sample(n = 300, seed = 1)

# OLS instantaneous structure without pruning (no extra packages required)
model <- lingam_varma(sample$data,
  order = c(1, 1), criterion = NULL,
  reg_method = "ols", prune = FALSE
)
model$causal_order
#> [1] 1 2 3
round(model$adjacency_matrices$psis[1, , ], 2) # instantaneous B0
#>      x0    x1 x2
#> x0 0.00  0.00  0
#> x1 0.61  0.00  0
#> x2 0.11 -0.56  0
```
