# Generate sample data with a latent confounder (for BottomUpParceLiNGAM)

Generates the 7-variable model used in the ParceLiNGAM tutorial, where
`x6` is an unobserved (latent) common cause of `x2` and `x3`. Only
`x0`-`x5` are returned as observed data; `x6` is not included.

## Usage

``` r
generate_parce_sample(n = 1000L, seed = NULL)
```

## Arguments

- n:

  number of samples (default: 1000)

- seed:

  random seed (default: NULL, i.e. do not reset the RNG state)

## Value

list with three elements:

- `data`: data.frame of the 6 observed variables (`x0`-`x5`).

- `adjacency_matrix`: the true 6x6 adjacency matrix among the observed
  variables, following the `m[to, from]` convention. The `x2`-`x3`
  entries (which share the latent confounder `x6` and have no direct
  edge between them) are `NA`, matching the convention used by
  [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md).

- `confounded_pair`: 1-based column positions of `x2` and `x3` (the pair
  sharing the latent confounder).

## Details

The data-generating process (all error terms are `Uniform(0, 1)`):


    x6 (latent) ~ Uniform(0, 1)
    x3 = 2.0 * x6 + e
    x2 = 2.0 * x6 + e
    x0 = 0.5 * x3 + e
    x1 = 0.5 * x0 + 0.5 * x2 + e
    x5 = 0.5 * x0 + e
    x4 = 0.5 * x0 - 0.5 * x2 + e

## Examples

``` r
confounded <- generate_parce_sample(n = 500, seed = 42)
head(confounded$data)
#>          x0       x1       x2       x3          x4        x5
#> 1 1.0369696 2.847403 2.677905 1.966117 -0.54667538 0.5788658
#> 2 1.9807394 2.397312 1.936897 2.051287  0.96611782 1.9233741
#> 3 0.5715207 1.681726 1.392124 1.091840  0.03568159 0.6347020
#> 4 2.1567712 3.067590 2.200256 2.472016  0.52004497 1.4901835
#> 5 1.0660913 2.258461 1.782511 1.398853 -0.19645547 1.4941921
#> 6 1.6641994 2.096731 1.060419 1.931614  0.99494363 1.3631726
confounded$adjacency_matrix
#>     x0 x1   x2  x3 x4 x5
#> x0 0.0  0  0.0 0.5  0  0
#> x1 0.5  0  0.5 0.0  0  0
#> x2 0.0  0  0.0  NA  0  0
#> x3 0.0  0   NA 0.0  0  0
#> x4 0.5  0 -0.5 0.0  0  0
#> x5 0.5  0  0.0 0.0  0  0
confounded$confounded_pair
#> [1] 3 4
```
