# ブートストラップ結果から因果効果の代表値の隣接行列を作成

ブートストラップ結果から因果効果の代表値の隣接行列を作成

## Usage

``` r
get_adjacency_matrix_summary(
  result,
  stat = "median",
  min_causal_effect = NULL,
  min_probability = NULL,
  labels = NULL
)
```

## Arguments

- result:

  BootstrapResult オブジェクト

- stat:

  代表値 ("mean" or "median")

- min_causal_effect:

  因果効果の最小閾値（これ以下はゼロ扱い）(NULL = 0)

- min_probability:

  この確率未満のエッジはゼロにする (NULL = 0)

- labels:

  変数名ベクトル (NULL可)

## Value

隣接行列 (n_features x n_features)。 **規則: `B[i, j]` は変数 j から変数
i への因果係数（j → i）。**
[`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
の `adjacency_matrix` と同じ規則。

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

bs_model <- lingam_direct_bootstrap(LiNGAM_sample_1000$data, n_sampling = 30L, seed = 42)
#> Bootstrap: 30 iterations, method=adaptive_lasso (sequential)
#>   iteration 1 / 30
#>   iteration 10 / 30
#>   iteration 20 / 30
#>   iteration 30 / 30
#> Completed in 2.1 seconds.
get_adjacency_matrix_summary(bs_model)
#>          [,1]       [,2]       [,3]     [,4]       [,5] [,6]
#> [1,] 0.000000 0.05299398  0.0000000 3.032919  0.1045919    0
#> [2,] 2.986914 0.00000000  2.0027022 0.000000  0.0000000    0
#> [3,] 0.000000 0.40422428  0.0000000 6.003637 -0.1387932    0
#> [4,] 0.000000 0.00000000  0.1616537 0.000000  0.0000000    0
#> [5,] 8.030490 0.90679690 -1.0165052 0.000000  0.0000000    0
#> [6,] 4.018861 0.00000000  0.0000000 0.000000  0.0000000    0
```
