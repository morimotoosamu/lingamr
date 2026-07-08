# Generate sample data for Multi-Group Direct LiNGAM (2 groups, 6 variables)

Generates two datasets that share the same causal structure as
[`generate_lingam_sample_6()`](https://morimotoosamu.github.io/lingamr/reference/generate_lingam_sample_6.md)
(`x3 -> x0, x3 -> x2, x0 -> x1, x2 -> x1, x0 -> x4, x2 -> x4, x0 -> x5`)
but with different structural coefficients per group, following the
multi-dataset tutorial's setup.

## Usage

``` r
generate_multi_group_sample(n = c(1000, 1000), seed = 42L)
```

## Arguments

- n:

  Numeric vector of length 2: sample size per group (default
  `c(1000, 1000)`).

- seed:

  Random seed (default 42). Group 2 uses an internally offset seed so
  the two groups are independently drawn.

## Value

A list with:

- `data_list`: named list (`group1`, `group2`) of data frames, each with
  columns `x0`..`x5`.

- `adjacency_matrices`: named list of the true adjacency matrix per
  group (`m[to, from] = coef`, same convention as
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)).

- `causal_order`: the true causal order shared by both groups (1-based
  indices into `x0`..`x5`).

## Examples

``` r
mg <- generate_multi_group_sample()
lapply(mg$data_list, head)
#> $group1
#>         x0        x1       x2        x3        x4        x5
#> 1 2.814924 18.017120 4.543655 0.6333728 18.160090 12.236660
#> 2 1.889685 10.956005 2.188091 0.3175366 13.172754  7.932657
#> 3 1.008905  6.990652 1.953131 0.2409218  6.702107  4.797122
#> 4 1.965690 12.296763 2.847148 0.3784141 13.224002  8.685252
#> 5 1.698178  9.698147 2.145058 0.3521443 11.673495  7.366258
#> 6 1.412372  8.640107 1.929980 0.2977585 10.024075  6.340899
#> 
#> $group2
#>          x0        x1        x2         x3        x4        x5
#> 1 0.7259014  5.482225 0.9301592 0.01259095  5.275903  3.321061
#> 2 2.4321051 17.252303 3.3989705 0.41696287 16.459975 11.368701
#> 3 1.5550457 10.342355 1.8713591 0.24518297 10.520165  7.859683
#> 4 3.7253048 26.770305 5.3418037 0.78548966 23.759223 16.826226
#> 5 0.9229987  4.845742 0.6255814 0.09606444  7.267970  4.910533
#> 6 2.3835023 18.944552 3.8564374 0.57442996 14.830625 11.001024
#> 
mg$adjacency_matrices$group1
#>    x0 x1 x2 x3 x4 x5
#> x0  0  0  0  3  0  0
#> x1  3  0  2  0  0  0
#> x2  0  0  0  6  0  0
#> x3  0  0  0  0  0  0
#> x4  8  0 -1  0  0  0
#> x5  4  0  0  0  0  0
```
