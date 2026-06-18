# Enumerate bootstrap paths between two variables in a VAR-LiNGAM model

Builds the time-expanded graph for every bootstrap sample and enumerates
all directed paths from the source (at `from_lag`) to the destination
(at `to_lag`), reporting each path's bootstrap probability and median
effect. Port of the Python reference `VARBootstrapResult.get_paths`.

## Usage

``` r
get_var_paths(
  result,
  from_index,
  to_index,
  from_lag = 0,
  to_lag = 0,
  min_causal_effect = NULL
)
```

## Arguments

- result:

  a VARBootstrapResult object

- from_index:

  source variable (1-based)

- to_index:

  destination variable (1-based)

- from_lag:

  lag of the source (default 0)

- to_lag:

  lag of the destination (default 0); must satisfy `to_lag <= from_lag`

- min_causal_effect:

  minimum \|effect\| threshold (NULL = 0)

## Value

a data frame (path, effect, probability), one row per distinct path

## Details

Node indices in the returned `path` are 1-based positions in the
time-expanded graph: column j of block L (lag L) corresponds to index
`n_features * L + j`.

## Examples

``` r
s <- generate_varlingam_sample(n = 500, seed = 42)
bs <- lingam_var_bootstrap(s$data,
  n_sampling = 10L, criterion = NULL,
  reg_method = "ols", prune = FALSE, seed = 1, verbose = FALSE
)
get_var_paths(bs, from_index = 1, to_index = 3)
#>      path        effect probability
#> 1 1, 2, 3 -0.3077611427           1
#> 2    1, 3 -0.0007678152           1
```
