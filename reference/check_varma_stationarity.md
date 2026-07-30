# Check the stationarity and invertibility of a fitted VARMA-LiNGAM model

Inspects the eigenvalues of the companion matrices of the reduced-form
AR coefficients (Phi, stationarity) and MA coefficients (Theta,
invertibility) stored in the result. The process is stationary when
every AR eigenvalue lies strictly inside the unit circle, and invertible
when every MA eigenvalue does; a modulus on or outside the circle
signals a (near-)unit root or a non-invertible MA polynomial, under
which the VARMA-LiNGAM estimates (and the residual filtering) are
unreliable. Invertibility is worth checking here because the
Hannan-Rissanen estimator does not enforce it.

## Usage

``` r
check_varma_stationarity(result, tol = 1)
```

## Arguments

- result:

  a `VARMALiNGAMResult` from
  [`lingam_varma()`](https://morimotoosamu.github.io/lingamr/reference/lingam_varma.md)

- tol:

  threshold for the eigenvalue moduli (default 1)

## Value

a `varma_stationarity` object (list) with `ar_moduli` / `ma_moduli`
(sorted descending; empty when p = 0 / q = 0), `max_ar_modulus`,
`max_ma_modulus`, `is_stationary`, `is_invertible`, `order`, and `tol`.

## References

Stationarity diagnostics in the spirit of the VARLiNGAM R code of
Moneta, A., Entner, D., Hoyer, P. O., & Coad, A. (2013), *Oxford
Bulletin of Economics and Statistics*, 75(5), 705-730.
<https://sites.google.com/site/dorisentner/publications/VARLiNGAM>

## Examples

``` r
s <- generate_varmalingam_sample(n = 1000, seed = 42)
m <- lingam_varma(s$data,
  order = c(1, 1), criterion = NULL,
  reg_method = "ols", prune = FALSE
)
check_varma_stationarity(m)
#> === VARMA Stationarity / Invertibility Check ===
#> Order (p, q):         (1, 1)
#> Max |AR eigenvalue|:  0.4763  (threshold 1.00)
#> Stationary:           YES
#> Max |MA eigenvalue|:  0.2974  (threshold 1.00)
#> Invertible:           YES
```
