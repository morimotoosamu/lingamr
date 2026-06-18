# Create an adjacency matrix of representative causal-effect values from bootstrap results

Create an adjacency matrix of representative causal-effect values from
bootstrap results

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

  BootstrapResult object

- stat:

  Representative statistic ("mean" or "median")

- min_causal_effect:

  Minimum threshold for the causal effect (values at or below this are
  treated as zero) (NULL = 0)

- min_probability:

  Edges below this probability are set to zero (NULL = 0)

- labels:

  Vector of variable names (NULL allowed)

## Value

Adjacency matrix (n_features x n_features). **Rule: `B[i, j]` is the
causal coefficient from variable j to variable i (j -\> i).** Same rule
as the `adjacency_matrix` of
[`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md).

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

bs_model <- lingam_direct_bootstrap(LiNGAM_sample_1000$data, n_sampling = 30L, seed = 42)
#> Bootstrap: 30 iterations, method=adaptive_lasso (sequential)
#>   iteration 1 / 30
#>   iteration 10 / 30
#>   iteration 20 / 30
#>   iteration 30 / 30
#> Completed in 1.1 seconds.
get_adjacency_matrix_summary(bs_model)
#>          [,1]       [,2]       [,3]     [,4]       [,5] [,6]
#> [1,] 0.000000 0.05299398  0.0000000 3.032919  0.1045919    0
#> [2,] 2.986914 0.00000000  2.0027022 0.000000  0.0000000    0
#> [3,] 0.000000 0.40422428  0.0000000 6.003637 -0.1387932    0
#> [4,] 0.000000 0.00000000  0.1616537 0.000000  0.0000000    0
#> [5,] 8.030490 0.90679690 -1.0165052 0.000000  0.0000000    0
#> [6,] 4.018861 0.00000000  0.0000000 0.000000  0.0000000    0
```
