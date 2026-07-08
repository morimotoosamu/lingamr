# Generate sample data with a latent confounder (for RCD)

Generates the 7-variable model used in the RCD tutorial, where `x6` is
an unobserved (latent) common cause of `x2` and `x4`. Only `x0`-`x5` are
returned as observed data; `x6` is not included.

## Usage

``` r
generate_rcd_sample(n = 300L, seed = NULL)
```

## Arguments

- n:

  number of samples (default: 300)

- seed:

  random seed (default: NULL, i.e. do not reset the RNG state)

## Value

list with four elements:

- `data`: data.frame of the 6 observed variables (`x0`-`x5`).

- `adjacency_matrix`: the true 6x6 adjacency matrix among the observed
  variables, following the `m[to, from]` convention. The `x2`-`x4`
  entries (which share the latent confounder `x6` and have no direct
  edge between them) are `NA`, matching the convention used by
  [`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md).

- `ancestors_list`: the true ancestor sets (1-based column positions),
  usable as a test oracle for
  [`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)'s
  `ancestors_list`.

- `confounded_pair`: 1-based column positions of `x2` and `x4` (the pair
  sharing the latent confounder).

## Details

The data-generating process (all error terms `e()` are super-Gaussian,
`rnorm(n, 0, 0.5)^3`):


    x5 ~ e();  x6 (latent) ~ e()
    x1 = 0.6 * x5 + e()
    x3 = 0.5 * x5 + e()
    x0 = 1.0 * x1 + 1.0 * x3 + e()
    x2 = 0.8 * x0 - 0.6 * x6 + e()
    x4 = 1.0 * x0 - 0.5 * x6 + e()

## Examples

``` r
confounded <- generate_rcd_sample(n = 300, seed = 1)
head(confounded$data)
#>            x0           x1         x2           x3         x4            x5
#> 1 -0.96839700 -0.023398020 -0.7514703 -0.473147026 -0.9674898 -0.0307310344
#> 2  1.31480945  0.424388510  1.0389688  0.001304296  1.4111680  0.0007741684
#> 3 -0.85973904 -0.025330472 -1.1731736 -0.034157640 -1.3389435 -0.0729373397
#> 4 -0.76463976  0.324413350 -0.7116352  0.078719765 -0.3333176  0.5074829194
#> 5  0.08152383  0.002364107 -0.2546673  0.036715357 -0.2040856  0.0044720536
#> 6 -0.30389184 -0.225029302 -0.5032249 -0.172267536  0.5221642 -0.0690391705
confounded$adjacency_matrix
#>     x0 x1 x2 x3 x4  x5
#> x0 0.0  1  0  1  0 0.0
#> x1 0.0  0  0  0  0 0.6
#> x2 0.8  0  0  0 NA 0.0
#> x3 0.0  0  0  0  0 0.0
#> x4 1.0  0 NA  0  0 0.0
#> x5 0.0  0  0  0  0 0.0
confounded$ancestors_list
#> $x0
#> [1] 2 4 6
#> 
#> $x1
#> [1] 6
#> 
#> $x2
#> [1] 1 2 4 6
#> 
#> $x3
#> [1] 6
#> 
#> $x4
#> [1] 1 2 4 6
#> 
#> $x5
#> integer(0)
#> 
confounded$confounded_pair
#> [1] 3 5
```
