# Fit a VARMA(p, q) model by two-stage Hannan-Rissanen

Stage 1 fits a long VAR(h) with intercept whose residuals estimate the
innovations; stage 2 regresses X on its own lags 1..p and the stage-1
residual lags 1..q (with intercept) by OLS.

## Usage

``` r
fit_varma_hr(X, p_order, q_order, h = NULL)
```

## Arguments

- X:

  numeric matrix (n_samples x n_features), rows ordered in time

- p_order:

  AR order (non-negative integer)

- q_order:

  MA order (non-negative integer)

- h:

  long AR order for stage 1 (NULL = use
  [`hr_long_ar_order()`](https://morimotoosamu.github.io/lingamr/reference/hr_long_ar_order.md);
  0 when q = 0, where stage 1 is skipped)

## Value

list with `phis` (array (p, m, m)), `thetas` (array (q, m, m)), `const`
(numeric m), `residuals` (stage-2 OLS residuals, one row per observation
in the estimation window t = h + q + 1 .. n), and `h`
