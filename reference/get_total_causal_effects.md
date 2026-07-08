# Get a list of total causal effects

Get a list of total causal effects

## Usage

``` r
get_total_causal_effects(result, min_causal_effect = NULL)
```

## Arguments

- result:

  BootstrapResult object

- min_causal_effect:

  Minimum threshold for the causal effect (NULL = 0)

## Value

data.frame (from, to, effect, probability)

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
#> Completed in 0.2 seconds.

get_total_causal_effects(bs_model)
#>    from to       effect probability
#> 1     1  6  4.021097096  1.00000000
#> 2     1  2  2.984272432  0.96666667
#> 3     1  5  8.009251257  0.96666667
#> 4     3  2  1.981440607  0.96666667
#> 5     3  5 -1.099920161  0.96666667
#> 6     4  1  3.032918694  0.96666667
#> 7     4  2 21.057962053  0.96666667
#> 8     4  3  6.003637462  0.96666667
#> 9     4  5 18.277681670  0.96666667
#> 10    4  6 12.184927850  0.96666667
#> 11    6  2  0.098273438  0.63333333
#> 12    5  2 -0.049991036  0.60000000
#> 13    3  6 -0.065251117  0.56666667
#> 14    6  5 -0.048877851  0.56666667
#> 15    3  1 -0.035631751  0.53333333
#> 16    1  3 -0.017162707  0.46666667
#> 17    5  6  0.033385392  0.43333333
#> 18    6  3  0.067423011  0.43333333
#> 19    2  5 -0.047501331  0.40000000
#> 20    2  6 -0.016891773  0.36666667
#> 21    1  4 -0.006219983  0.03333333
#> 22    2  1  0.147945030  0.03333333
#> 23    2  3  0.278509201  0.03333333
#> 24    2  4  0.046110069  0.03333333
#> 25    3  4  0.135177264  0.03333333
#> 26    5  1  0.103731652  0.03333333
#> 27    5  3 -0.142910453  0.03333333
#> 28    5  4 -0.010957236  0.03333333
#> 29    6  4 -0.004996343  0.03333333
```
