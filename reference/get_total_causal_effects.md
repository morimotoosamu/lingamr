# 総合因果効果リストを取得

総合因果効果リストを取得

## Usage

``` r
get_total_causal_effects(result, min_causal_effect = NULL)
```

## Arguments

- result:

  BootstrapResult オブジェクト

- min_causal_effect:

  因果効果の最小閾値 (NULL = 0)

## Value

data.frame (from, to, effect, probability)

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

bs_model <- lingam_direct_bootstrap(LiNGAM_sample_1000$data, n_sampling = 30L, seed = 42)
#> Bootstrap: 30 iterations, method=adaptive_lasso (sequential)
#>   iteration 1 / 30
#>   iteration 10 / 30
#>   iteration 20 / 30
#>   iteration 30 / 30
#> Completed in 2.0 seconds.

get_total_causal_effects(bs_model)
#>    from to      effect probability
#> 1     1  6  4.01907390  1.00000000
#> 2     1  2  2.93210481  0.96666667
#> 3     1  5  8.00575712  0.96666667
#> 4     3  2  1.94023398  0.96666667
#> 5     3  5 -1.18014784  0.96666667
#> 6     4  1  3.03291869  0.96666667
#> 7     4  2 21.05796205  0.96666667
#> 8     4  3  6.00363746  0.96666667
#> 9     4  5 18.27768167  0.96666667
#> 10    4  6 12.18492785  0.96666667
#> 11    6  2  0.19623886  0.06666667
#> 12    2  1  0.14794503  0.03333333
#> 13    2  3  0.27850920  0.03333333
#> 14    2  4  0.04611007  0.03333333
#> 15    2  5  0.90679690  0.03333333
#> 16    2  6  0.59359217  0.03333333
#> 17    3  4  0.16191585  0.03333333
#> 18    3  6 -0.35827083  0.03333333
#> 19    5  1  0.10506715  0.03333333
#> 20    5  3 -0.13869103  0.03333333
#> 21    5  6  0.41846402  0.03333333
```
