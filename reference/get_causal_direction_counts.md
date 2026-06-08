# 因果方向のカウント・割合・因果効果を取得

因果方向のカウント・割合・因果効果を取得

## Usage

``` r
get_causal_direction_counts(
  result,
  n_directions = NULL,
  min_causal_effect = NULL,
  split_by_causal_effect_sign = FALSE,
  labels = NULL
)
```

## Arguments

- result:

  BootstrapResult オブジェクト

- n_directions:

  上位何件を返すか (NULL = 全て)

- min_causal_effect:

  因果効果の最小閾値 (NULL = 0)

- split_by_causal_effect_sign:

  因果効果の符号で分割するか

- labels:

  変数名ベクトル (NULL可。指定するとfrom_name, to_name列を追加)

## Value

A data frame containing the following columns:

- `from`, `to`: 1-based indices of the causal (from) and effect (to)
  variables.

- `count`: Number of bootstrap samples in which this specific causal
  direction was identified.

- `proportion`: The frequency of the direction's occurrence (count /
  n_sampling), representing its bootstrap probability.

- `mean_effect`: The average value of the estimated causal effects
  across samples where this direction was identified.

- `median_effect`: The median value of the estimated causal effects,
  providing a robust estimate of the effect size.

- `sd_effect`: The standard deviation of the causal effect estimates,
  indicating the stability of the effect size.

- `ci_lower`, `ci_upper`: The lower (2.5%) and upper (97.5%) bounds of
  the bootstrap confidence interval for the causal effect.

- `sign` (optional): The sign of the causal effect (1 for positive, -1
  for negative), included if `split_by_causal_effect_sign = TRUE`.

- `from_name`, `to_name` (optional): Character labels for the variables,
  included if `labels` were provided.

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

bs_model <- lingam_direct_bootstrap(LiNGAM_sample_1000$data, n_sampling = 30L, seed = 42)
#> Bootstrap: 30 iterations, method=adaptive_lasso (sequential)
#>   iteration 1 / 30
#>   iteration 10 / 30
#>   iteration 20 / 30
#>   iteration 30 / 30
#> Completed in 1.8 seconds.

get_causal_direction_counts(bs_model, labels = names(LiNGAM_sample_1000$data))
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
#>       ci_upper from_name to_name
#> 1   4.03756004        x0      x5
#> 2   3.04468103        x0      x1
#> 3   8.08567538        x0      x4
#> 4   2.02581141        x2      x1
#> 5  -0.98775234        x2      x4
#> 6   3.09090167        x3      x0
#> 7   6.05479799        x3      x2
#> 8   0.05299398        x1      x0
#> 9   0.40422428        x1      x2
#> 10  0.90679690        x1      x4
#> 11  0.16165370        x2      x3
#> 12  0.10459193        x4      x0
#> 13 -0.13879324        x4      x2
```
