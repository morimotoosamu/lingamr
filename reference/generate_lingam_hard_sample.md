# Generate a challenging sample data for Direct LiNGAM

Generates a dataset with conditions that make causal estimation
difficult:

1.  High multicollinearity among predictors

2.  Moderate sample size relative to variables

3.  True coefficients of similar magnitude

## Usage

``` r
generate_lingam_hard_sample(n = 200L, seed = 42L, collinearity = 0.95)
```

## Arguments

- n:

  number of samples (default: 200)

- seed:

  random seed (default: 42)

- collinearity:

  strength of multicollinearity (0 to 1, default: 0.95)

## Value

list(data, true_adjacency)

## Details

These conditions destabilize OLS initial estimates in Adaptive LASSO,
making Ridge-initialized Adaptive LASSO preferable.

## Examples

``` r
result <- generate_lingam_hard_sample()
result$true_adjacency
#>     x0  x1  x2 x3 x4 x5 x6 x7 x8
#> x0 0.0 0.0 0.0  0  0  0  0  0  0
#> x1 0.0 0.0 0.0  0  0  0  0  0  0
#> x2 0.0 0.0 0.0  0  0  0  0  0  0
#> x3 0.0 0.0 0.0  0  0  0  0  0  0
#> x4 0.0 0.0 0.0  0  0  0  0  0  0
#> x5 1.5 1.5 1.5  0  0  0  0  0  0
#> x6 0.0 1.0 1.0  1  1  0  0  0  0
#> x7 2.0 0.0 0.0  2  0  0  0  0  0
#> x8 0.0 0.0 0.0  0  0  1  1  0  0
head(result$data)
#>          x0        x1        x2        x3        x4       x5       x6       x7
#> 1 1.0228248 0.8886297 0.9693088 0.9349473 0.7947395 5.298107 4.166792 4.281243
#> 2 0.4830486 0.4748335 0.2787734 0.2895976 0.2664368 2.228899 1.343471 2.033147
#> 3 0.3727160 0.3013776 0.4418671 0.3585968 0.4657222 2.435443 1.570280 1.732436
#> 4 1.1951289 1.1560823 1.1158846 1.0539815 1.0437446 6.023135 4.747024 5.105043
#> 5 0.4059819 0.3034858 0.2156491 0.3155540 0.2783906 1.961220 1.356371 1.488881
#> 6 0.9547901 0.9622644 0.8374879 0.8856775 0.9972516 4.823226 3.801561 3.977539
#>          x8
#> 1 10.173626
#> 2  4.010030
#> 3  4.205728
#> 4 11.537226
#> 5  3.830753
#> 6  8.669491
```
