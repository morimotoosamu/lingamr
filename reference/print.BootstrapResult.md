# BootstrapResult の内容を表示

BootstrapResult の内容を表示

## Usage

``` r
# S3 method for class 'BootstrapResult'
print(x, ...)
```

## Arguments

- x:

  BootstrapResult オブジェクト

- ...:

  追加の引数 (S3メソッド互換用)

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

bs_model <- LiNGAM_sample_1000$data |>
  lingam_direct_bootstrap(n_sampling = 30L, seed = 42)
#> Bootstrap: 30 iterations, method=adaptive_lasso (sequential)
#>   iteration 1 / 30
#>   iteration 10 / 30
#>   iteration 20 / 30
#>   iteration 30 / 30
#> Completed in 1.9 seconds.

bs_model |>
  print()
#> BootstrapResult: 30 samplings, 6 features
```
