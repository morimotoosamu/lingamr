# Residual matrix to diagnose for a VAR-LiNGAM model

Returns the series targeted by the residual diagnostics: either the
LiNGAM innovations `e_t = (I - B0) n_t` (the independent errors) or the
reduced-form VAR residuals `n_t`. Shared by the normality tests and the
QQ plot.

## Usage

``` r
compute_varlingam_residuals(result, on = c("innovations", "var"))
```

## Arguments

- result:

  a `VARLiNGAMResult`

- on:

  "innovations" or "var"

## Value

residual matrix (n_obs x n_features), column names preserved
