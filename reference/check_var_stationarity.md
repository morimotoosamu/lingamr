# Check the stationarity of a fitted VAR-LiNGAM model

Recovers the reduced-form VAR coefficients `M_k = (I - B0)^{-1} B_k`
from the structural matrices and inspects the eigenvalues of the VAR
companion matrix. The process is stationary when every eigenvalue lies
strictly inside the unit circle (all moduli \< 1); a modulus on or
outside it signals a (near-)unit root, under which the VAR-LiNGAM
estimates are unreliable.

## Usage

``` r
check_var_stationarity(result, tol = 1)
```

## Arguments

- result:

  a `VARLiNGAMResult` from
  [`lingam_var()`](https://morimotoosamu.github.io/lingamr/reference/lingam_var.md)

- tol:

  stationarity threshold for the eigenvalue moduli (default 1)

## Value

a `var_stationarity` object (list) with `moduli` (sorted descending),
`max_modulus`, `is_stationary` (logical), `lags`, and `tol`.

## References

Stationarity diagnostics in the spirit of the VARLiNGAM R code of
Moneta, A., Entner, D., Hoyer, P. O., & Coad, A. (2013), *Oxford
Bulletin of Economics and Statistics*, 75(5), 705-730.
<https://sites.google.com/site/dorisentner/publications/VARLiNGAM>

## Examples

``` r
s <- generate_varlingam_sample(n = 1000, seed = 42)
m <- lingam_var(s$data, lags = 1, reg_method = "ols", prune = FALSE)
check_var_stationarity(m)
#> === VAR Stationarity Check ===
#> Lag order:         1
#> Max |eigenvalue|:  0.5038  (threshold 1.00)
#> Stationary:        YES
```
