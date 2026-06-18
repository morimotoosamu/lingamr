# Get bootstrap probabilities

Get bootstrap probabilities

## Usage

``` r
get_probabilities(result, min_causal_effect = NULL)
```

## Arguments

- result:

  BootstrapResult object

- min_causal_effect:

  Minimum threshold for the causal effect (NULL = 0)

## Value

Probability matrix (n_features x n_features)

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

bs_model <- lingam_direct_bootstrap(LiNGAM_sample_1000$data, n_sampling = 30L, seed = 42)
#> Bootstrap: 30 iterations, method=adaptive_lasso (sequential)
#>   iteration 1 / 30
#>   iteration 10 / 30
#>   iteration 20 / 30
#>   iteration 30 / 30
#> Completed in 0.7 seconds.

get_probabilities(bs_model)
#>           [,1]       [,2]       [,3]      [,4]       [,5] [,6]
#> [1,] 0.0000000 0.03333333 0.00000000 0.9666667 0.03333333    0
#> [2,] 0.9666667 0.00000000 0.96666667 0.0000000 0.00000000    0
#> [3,] 0.0000000 0.03333333 0.00000000 0.9666667 0.03333333    0
#> [4,] 0.0000000 0.00000000 0.03333333 0.0000000 0.00000000    0
#> [5,] 0.9666667 0.03333333 0.96666667 0.0000000 0.00000000    0
#> [6,] 1.0000000 0.00000000 0.00000000 0.0000000 0.00000000    0
```
