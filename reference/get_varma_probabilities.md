# Bootstrap probabilities for a VARMA-LiNGAM model

Returns, for each entry of the joined matrix, the fraction of bootstrap
samples in which that edge exceeded `min_causal_effect`.

## Usage

``` r
get_varma_probabilities(result, min_causal_effect = NULL)
```

## Arguments

- result:

  a VARMABootstrapResult object

- min_causal_effect:

  minimum \|effect\| threshold (NULL = 0)

## Value

probability matrix (n_features x n_features\*(1 + p + q)). Columns
1..n_features are the instantaneous block, the next p blocks are the AR
lags 1..p (psi), and the final q blocks are the MA terms 1..q (omega).
`P[i, j]` is the probability of the edge j -\> i.

## Examples

``` r
s <- generate_varmalingam_sample(n = 300, seed = 42)
bs <- lingam_varma_bootstrap(s$data,
  n_sampling = 5L, order = c(1, 1), criterion = NULL,
  reg_method = "ols", prune = FALSE, seed = 1, verbose = FALSE
)
get_varma_probabilities(bs)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#> [1,]    0    0    0    1    1    1    1    1    1
#> [2,]    1    0    0    1    1    1    1    1    1
#> [3,]    1    1    0    1    1    1    1    1    1
```
