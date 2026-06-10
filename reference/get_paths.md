# 指定した2変数間の全パスとブートストラップ確率を取得

指定した2変数間の全パスとブートストラップ確率を取得

## Usage

``` r
get_paths(result, from_index, to_index, min_causal_effect = NULL)
```

## Arguments

- result:

  BootstrapResult オブジェクト

- from_index:

  始点インデックス (1-based)

- to_index:

  終点インデックス (1-based)

- min_causal_effect:

  因果効果の最小閾値 (NULL = 0)

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
#> Completed in 2.0 seconds.
get_paths(bs_model, 1, 6)
#>   path   effect probability
#> 1 1, 6 4.018861           1
```
