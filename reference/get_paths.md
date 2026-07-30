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

bs_model <- lingam_direct_bootstrap(LiNGAM_sample_1000$data,
  n_sampling = 30L, reg_method = "ols", seed = 42
)
#> Bootstrap: 30 iterations, method=ols (sequential)
#>   iteration 1 / 30
#>   iteration 10 / 30
#>   iteration 20 / 30
#>   iteration 30 / 30
#> Completed in 0.1 seconds.
get_paths(bs_model, 1, 6)
#>           path        effect probability
#> 1         1, 6  4.0038996964  1.00000000
#> 2      1, 5, 6  0.2347844496  0.40000000
#> 3      1, 2, 6 -0.0533555108  0.33333333
#> 4   1, 2, 5, 6 -0.0027424587  0.16666667
#> 5   1, 5, 2, 6  0.0052283471  0.13333333
#> 6 1, 3, 2,....  0.0002273752  0.03333333
#> 7   1, 3, 2, 6  0.0016654412  0.03333333
#> 8   1, 3, 5, 6  0.0018510027  0.03333333
#> 9      1, 3, 6 -0.0051552729  0.03333333
```
