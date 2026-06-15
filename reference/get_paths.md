# Get all paths between two specified variables and their bootstrap probabilities

Get all paths between two specified variables and their bootstrap
probabilities

## Usage

``` r
get_paths(result, from_index, to_index, min_causal_effect = NULL)
```

## Arguments

- result:

  BootstrapResult object

- from_index:

  Start index (1-based)

- to_index:

  End index (1-based)

- min_causal_effect:

  Minimum threshold for the causal effect (NULL = 0)

## Value

data.frame (path, effect, probability)

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

bs_model <- lingam_direct_bootstrap(LiNGAM_sample_1000$data, n_sampling = 30L, seed = 42)
#> Bootstrap: 30 iterations, method=adaptive_lasso (sequential)
#>   iteration 1 / 30
#>   iteration 10 / 30
#>   iteration 20 / 30
#>   iteration 30 / 30
#> Completed in 0.9 seconds.
get_paths(bs_model, 1, 6)
#>   path   effect probability
#> 1 1, 6 4.018861           1
```
