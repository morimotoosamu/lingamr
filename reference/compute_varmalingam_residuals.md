# Residual matrix to diagnose for a VARMA-LiNGAM model

Returns the series targeted by the residual diagnostics: either the
LiNGAM innovations `e_t = (I - B0) n_t` (the independent errors) or the
reduced-form VARMA residuals `n_t`. Shared by the normality tests and
the QQ plot.

## Usage

``` r
compute_varmalingam_residuals(result, on = c("innovations", "varma"))
```

## Arguments

- result:

  a `VARMALiNGAMResult`

- on:

  "innovations" or "varma"

## Value

residual matrix (n_obs x n_features), column names preserved
