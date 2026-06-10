# BootstrapResult を tidy な data.frame に変換

各因果方向の出現回数・割合・効果量の要約を返す。内部で
[`get_causal_direction_counts()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_direction_counts.md)
を呼び出すため、同関数の引数を `...` で渡せる。

## Usage

``` r
# S3 method for class 'BootstrapResult'
tidy(x, ...)
```

## Arguments

- x:

  [`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md)
  の返り値（`BootstrapResult` オブジェクト）

- ...:

  [`get_causal_direction_counts()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_direction_counts.md)
  に渡す引数 (`n_directions`, `min_causal_effect`,
  `split_by_causal_effect_sign`, `labels` など)

## Value

data.frame (from, to, count, proportion, ...)

## Examples

``` r
dat <- generate_lingam_sample_6()
bs <- lingam_direct_bootstrap(dat$data, n_sampling = 30L, seed = 42)
#> Bootstrap: 30 iterations, method=adaptive_lasso (sequential)
#>   iteration 1 / 30
#>   iteration 10 / 30
#>   iteration 20 / 30
#>   iteration 30 / 30
#> Completed in 2.0 seconds.
tidy(bs)
#>    from to count proportion mean_effect median_effect   sd_effect    ci_lower
#> 1     1  6    30 1.00000000  4.02010918    4.01886128 0.009510235  4.00298581
#> 2     1  2    29 0.96666667  2.98872289    2.98691378 0.029595289  2.94808621
#> 3     1  5    29 0.96666667  8.02805949    8.03049041 0.030493711  7.97909377
#> 4     3  2    29 0.96666667  2.00133847    2.00270222 0.015820048  1.97089548
#> 5     3  5    29 0.96666667 -1.01545079   -1.01650518 0.015564566 -1.03972182
#> 6     4  1    29 0.96666667  3.03159775    3.03291869 0.035791478  2.96698487
#> 7     4  3    29 0.96666667  6.00046795    6.00363746 0.031652501  5.94254816
#> 8     2  1     1 0.03333333  0.05299398    0.05299398 0.000000000  0.05299398
#> 9     2  3     1 0.03333333  0.40422428    0.40422428 0.000000000  0.40422428
#> 10    2  5     1 0.03333333  0.90679690    0.90679690 0.000000000  0.90679690
#> 11    3  4     1 0.03333333  0.16165370    0.16165370 0.000000000  0.16165370
#> 12    5  1     1 0.03333333  0.10459193    0.10459193 0.000000000  0.10459193
#> 13    5  3     1 0.03333333 -0.13879324   -0.13879324 0.000000000 -0.13879324
#>       ci_upper
#> 1   4.03756004
#> 2   3.04468103
#> 3   8.08567538
#> 4   2.02581141
#> 5  -0.98775234
#> 6   3.09090167
#> 7   6.05479799
#> 8   0.05299398
#> 9   0.40422428
#> 10  0.90679690
#> 11  0.16165370
#> 12  0.10459193
#> 13 -0.13879324
```
