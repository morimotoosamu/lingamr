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

bs_model <- lingam_direct_bootstrap(LiNGAM_sample_1000$data,
  n_sampling = 30L, reg_method = "ols", seed = 42
)
#> Bootstrap: 30 iterations, method=ols (sequential)
#>   iteration 1 / 30
#>   iteration 10 / 30
#>   iteration 20 / 30
#>   iteration 30 / 30
#> Completed in 0.4 seconds.
get_adjacency_matrix_summary(bs_model)
#>            [,1]         [,2]        [,3]         [,4]         [,5]         [,6]
#> [1,]  0.0000000  0.042882297 -0.03563175  3.091442871  0.107583403  0.000000000
#> [2,]  3.0682218  0.000000000  1.96645354 -0.008653801 -0.050117494  0.008629175
#> [3,] -0.2628120  0.408099958  0.00000000  6.006990357 -0.142910453  0.067423011
#> [4,]  0.0111159  0.001014947  0.13602291  0.000000000  0.009409755 -0.004996343
#> [5,]  8.0064511 -0.047501331 -1.05335717  0.382649748  0.000000000 -0.001975088
#> [6,]  4.0038997 -0.016891773  0.08530374 -0.079885472  0.029746764  0.000000000
```
