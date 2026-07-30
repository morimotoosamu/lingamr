# Generate sample data with unobserved variables (for CAM-UV)

Generates a 6-variable nonlinear additive model with two *unobserved*
variables, following the structure of the CAM-UV tutorial of the Python
`lingam` package: `u1` is an unobserved common cause of `x3` and `x4`
(an unobserved backdoor path, UBP), and `u2` is an unobserved
intermediate variable on the path from `x2` to `x5` (an unobserved
causal path, UCP). Only `x0`-`x5` are returned as observed data.

## Usage

``` r
generate_camuv_sample(n = 500L, seed = NULL)
```

## Arguments

- n:

  number of samples (default: 500)

- seed:

  random seed (default: NULL, i.e. do not reset the RNG state)

## Value

list with three elements:

- `data`: data.frame of the 6 observed variables (`x0`-`x5`).

- `adjacency_matrix`: the true 6x6 adjacency matrix among the observed
  variables, following the `m[to, from]` convention. Entries are 0/1
  edge indicators (the causal functions are nonlinear), directly
  comparable to the output of
  [`lingam_camuv()`](https://morimotoosamu.github.io/lingamr/reference/lingam_camuv.md).
  The `x3`-`x4` (UBP) and `x2`-`x5` (UCP) entries are `NA`, matching the
  convention used by
  [`lingam_camuv()`](https://morimotoosamu.github.io/lingamr/reference/lingam_camuv.md).

- `confounded_pairs`: 2-column integer matrix of the variable pairs
  (1-based column positions) connected through an unobserved variable,
  usable as a test oracle for
  [`lingam_camuv()`](https://morimotoosamu.github.io/lingamr/reference/lingam_camuv.md)'s
  `confounded_pairs`.

## Details

The data-generating process (all error terms `e()` are
`runif(n, -2.5, 2.5)`; the tutorial's per-call random constants are
replaced by the fixed constants below for reproducibility; each column
is standardized before being used as a cause, as in the tutorial):


    u1 (latent) ~ e();  x0..x5 ~ e()
    x3 = x3 + (u1 + 1.5)^2
    x4 = x4 + (u1 + 1.2)^2
    x1 = x1 + (x0 + 1.2)^2
    x3 = x3 + (x0 - 1.5)^2
    x4 = x4 + (x2 + 1.0)^2
    u2 (latent) = (x2 - 1.2)^2 + e()
    x5 = x5 + (u2 + 1.5)^2

## Examples

``` r
confounded <- generate_camuv_sample(n = 200, seed = 1)
head(confounded$data)
#>           x0         x1         x2       x3        x4         x5
#> 1 -0.8716163 -0.3465264  0.5446793 1.942204 1.2914382 0.30235032
#> 2 -0.4753220 -0.2896798 -1.0796731 2.666484 0.4709917 2.26086593
#> 3  0.2707997  0.7383932  1.5577423 1.990869 4.3486061 0.85648319
#> 4  1.5173293  2.0465938  1.3639420 1.307915 2.7211432 0.73327140
#> 5 -1.1088636 -0.5210769  1.5211244 3.259314 2.3822889 0.04943462
#> 6  1.4808349  2.3919707  0.7668779 2.224420 2.6211662 0.35048326
confounded$adjacency_matrix
#>    x0 x1 x2 x3 x4 x5
#> x0  0  0  0  0  0  0
#> x1  1  0  0  0  0  0
#> x2  0  0  0  0  0 NA
#> x3  1  0  0  0 NA  0
#> x4  0  0  1 NA  0  0
#> x5  0  0 NA  0  0  0
confounded$confounded_pairs
#>      var1 var2
#> [1,]    3    6
#> [2,]    4    5
```
