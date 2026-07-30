# Get DAG counts

Get DAG counts

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

  BootstrapResult object

- n_dags:

  How many of the top entries to return (NULL = all)

- min_causal_effect:

  Minimum threshold for the causal effect (NULL = 0)

- split_by_causal_effect_sign:

  Whether to split by the sign of the causal effect

## Value

list(dag = list of data.frames, count = integer vector)

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
#> Completed in 0.1 seconds.

get_directed_acyclic_graph_counts(bs_model)
#> $dag
#> $dag[[1]]
#>    from to
#> 1     1  2
#> 2     1  3
#> 3     1  5
#> 4     1  6
#> 5     3  2
#> 6     3  5
#> 7     4  1
#> 8     4  2
#> 9     4  3
#> 10    4  5
#> 11    4  6
#> 12    5  2
#> 13    6  2
#> 14    6  3
#> 15    6  5
#> 
#> $dag[[2]]
#>    from to
#> 1     1  2
#> 2     1  5
#> 3     1  6
#> 4     2  6
#> 5     3  1
#> 6     3  2
#> 7     3  5
#> 8     3  6
#> 9     4  1
#> 10    4  2
#> 11    4  3
#> 12    4  5
#> 13    4  6
#> 14    5  2
#> 15    5  6
#> 
#> $dag[[3]]
#>    from to
#> 1     1  2
#> 2     1  5
#> 3     1  6
#> 4     2  5
#> 5     2  6
#> 6     3  1
#> 7     3  2
#> 8     3  5
#> 9     3  6
#> 10    4  1
#> 11    4  2
#> 12    4  3
#> 13    4  5
#> 14    4  6
#> 15    5  6
#> 
#> $dag[[4]]
#>    from to
#> 1     1  2
#> 2     1  3
#> 3     1  5
#> 4     1  6
#> 5     2  5
#> 6     3  2
#> 7     3  5
#> 8     4  1
#> 9     4  2
#> 10    4  3
#> 11    4  5
#> 12    4  6
#> 13    6  2
#> 14    6  3
#> 15    6  5
#> 
#> $dag[[5]]
#>    from to
#> 1     1  2
#> 2     1  5
#> 3     1  6
#> 4     3  1
#> 5     3  2
#> 6     3  5
#> 7     3  6
#> 8     4  1
#> 9     4  2
#> 10    4  3
#> 11    4  5
#> 12    4  6
#> 13    5  2
#> 14    5  6
#> 15    6  2
#> 
#> $dag[[6]]
#>    from to
#> 1     1  2
#> 2     1  5
#> 3     1  6
#> 4     3  1
#> 5     3  2
#> 6     3  5
#> 7     3  6
#> 8     4  1
#> 9     4  2
#> 10    4  3
#> 11    4  5
#> 12    4  6
#> 13    5  2
#> 14    6  2
#> 15    6  5
#> 
#> $dag[[7]]
#>    from to
#> 1     1  4
#> 2     1  6
#> 3     2  1
#> 4     2  3
#> 5     2  4
#> 6     2  5
#> 7     2  6
#> 8     3  1
#> 9     3  4
#> 10    3  6
#> 11    5  1
#> 12    5  3
#> 13    5  4
#> 14    5  6
#> 15    6  4
#> 
#> $dag[[8]]
#>    from to
#> 1     1  2
#> 2     1  5
#> 3     1  6
#> 4     2  5
#> 5     3  1
#> 6     3  2
#> 7     3  5
#> 8     3  6
#> 9     4  1
#> 10    4  2
#> 11    4  3
#> 12    4  5
#> 13    4  6
#> 14    6  2
#> 15    6  5
#> 
#> $dag[[9]]
#>    from to
#> 1     1  2
#> 2     1  5
#> 3     1  6
#> 4     2  5
#> 5     2  6
#> 6     3  1
#> 7     3  2
#> 8     3  5
#> 9     3  6
#> 10    4  1
#> 11    4  2
#> 12    4  3
#> 13    4  5
#> 14    4  6
#> 15    6  5
#> 
#> $dag[[10]]
#>    from to
#> 1     1  2
#> 2     1  3
#> 3     1  5
#> 4     1  6
#> 5     2  5
#> 6     2  6
#> 7     3  2
#> 8     3  5
#> 9     3  6
#> 10    4  1
#> 11    4  2
#> 12    4  3
#> 13    4  5
#> 14    4  6
#> 15    5  6
#> 
#> 
#> $count
#>  [1] 9 4 4 4 3 2 1 1 1 1
#> 
```
