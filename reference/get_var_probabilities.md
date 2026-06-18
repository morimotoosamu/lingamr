# Bootstrap probabilities for a VAR-LiNGAM model

Returns, for each entry of the joined adjacency matrix, the fraction of
bootstrap samples in which that edge exceeded `min_causal_effect`.

## Usage

``` r
get_var_probabilities(result, min_causal_effect = NULL)
```

## Arguments

- result:

  a VARBootstrapResult object

- min_causal_effect:

  minimum \|effect\| threshold (NULL = 0)

## Value

probability matrix (n_features x n_features\*(1 + lags)). Columns
1..n_features are the instantaneous block; the next n_features are lag
1; etc. `P[i, j]` is the probability of the edge j -\> i.

## Examples

``` r
s <- generate_varlingam_sample(n = 500, seed = 42)
bs <- lingam_var_bootstrap(s$data,
  n_sampling = 10L, criterion = NULL,
  reg_method = "ols", prune = FALSE, seed = 1, verbose = FALSE
)
get_var_probabilities(bs)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    0    0    0    1    1    1
#> [2,]    1    0    0    1    1    1
#> [3,]    1    1    0    1    1    1
```
