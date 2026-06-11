# DAG カウントを取得

DAG カウントを取得

## Usage

``` r
get_directed_acyclic_graph_counts(
  result,
  n_dags = NULL,
  min_causal_effect = NULL,
  split_by_causal_effect_sign = FALSE
)
```

## Arguments

- result:

  BootstrapResult オブジェクト

- n_dags:

  上位何件を返すか (NULL = 全て)

- min_causal_effect:

  因果効果の最小閾値 (NULL = 0)

- split_by_causal_effect_sign:

  因果効果の符号で分割するか

## Value

list(dag = list of data.frames, count = integer vector)

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

bs_model <- lingam_direct_bootstrap(LiNGAM_sample_1000$data, n_sampling = 30L, seed = 42)
#> Bootstrap: 30 iterations, method=adaptive_lasso (sequential)
#>   iteration 1 / 30
#>   iteration 10 / 30
#>   iteration 20 / 30
#>   iteration 30 / 30
#> Completed in 1.4 seconds.

get_directed_acyclic_graph_counts(bs_model)
#> $dag
#> $dag[[1]]
#>   from to
#> 1    1  2
#> 2    1  5
#> 3    1  6
#> 4    3  2
#> 5    3  5
#> 6    4  1
#> 7    4  3
#> 
#> $dag[[2]]
#>   from to
#> 1    1  6
#> 2    2  1
#> 3    2  3
#> 4    2  5
#> 5    3  4
#> 6    5  1
#> 7    5  3
#> 
#> 
#> $count
#> [1] 29  1
#> 
```
