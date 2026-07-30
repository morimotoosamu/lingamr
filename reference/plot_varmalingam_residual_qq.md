# Q-Q plots of VARMA-LiNGAM residuals

Draws per-variable normal Q-Q plots of the residuals (analogous to the
Moneta `Gauss_Stats` visual check). Deviations from the reference line
indicate non-Gaussianity, which supports the LiNGAM assumption. Requires
ggplot2.

## Usage

``` r
plot_varmalingam_residual_qq(
  result,
  on = c("innovations", "varma"),
  ncol = 3,
  nrow = NULL
)
```

## Arguments

- result:

  a `VARMALiNGAMResult` from
  [`lingam_varma()`](https://morimotoosamu.github.io/lingamr/reference/lingam_varma.md)

- on:

  which series to plot: "innovations" (default) or "varma"

- ncol:

  number of facet columns

- nrow:

  number of facet rows (NULL = automatic)

## Value

a ggplot object

## References

Analogous to the residual visual check (Gauss_Stats) in the VARLiNGAM R
code of Moneta, A., Entner, D., Hoyer, P. O., & Coad, A. (2013), *Oxford
Bulletin of Economics and Statistics*, 75(5), 705-730.
<https://sites.google.com/site/dorisentner/publications/VARLiNGAM>

## Examples

``` r
s <- generate_varmalingam_sample(n = 1000, seed = 42)
m <- lingam_varma(s$data,
  order = c(1, 1), criterion = NULL,
  reg_method = "ols", prune = FALSE
)
# \donttest{
plot_varmalingam_residual_qq(m)

# }
```
