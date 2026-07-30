# Generate sample data from a VARMA-LiNGAM model

Generates a 3-variable time series following a VARMA(1,1)-LiNGAM model
with a strictly acyclic instantaneous structure B0, a lag-1
autoregressive matrix Phi1, a lag-1 moving-average matrix Theta1, and
non-Gaussian (uniform) errors. The reduced-form recursion is
`x(t) = Phi1 x(t-1) + n(t) + Theta1 n(t-1)` with
`n(t) = (I - B0)^{-1} e(t)`.

## Usage

``` r
generate_varmalingam_sample(n = 1000, seed = NULL)
```

## Arguments

- n:

  number of time points to return (after burn-in)

- seed:

  random seed (NULL allowed)

## Value

list with `data` (data frame, n x 3), the reduced-form true matrices
`true_B0`, `true_phi1`, `true_theta1`, and the structural-form
counterparts `true_psi1 = (I - B0) Phi1` and
`true_omega1 = (I - B0) Theta1 (I - B0)^{-1}` for comparison against
[`lingam_varma()`](https://morimotoosamu.github.io/lingamr/reference/lingam_varma.md)
estimates

## Examples

``` r
sample <- generate_varmalingam_sample(n = 500, seed = 1)
head(sample$data)
#>           x0         x1         x2
#> 1  0.1735288  1.2799008 -1.2155499
#> 2 -0.4072948  1.2166111 -1.4989248
#> 3 -0.9548545  0.5966513 -1.4501052
#> 4  0.2940208  1.0468406 -0.1676003
#> 5  0.5680147  1.1038743  0.7978025
#> 6 -0.1768722 -0.7770469  1.5696898
```
