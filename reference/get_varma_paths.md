# Enumerate bootstrap paths between two variables in a VARMA-LiNGAM model

Builds the time-expanded graph for every bootstrap sample and enumerates
all directed paths from the source (at `from_lag`) to the destination
(at `to_lag`), reporting each path's bootstrap probability and median
effect. Port of the Python reference `VARMABootstrapResult.get_paths`.

## Usage

``` r
get_varma_paths(
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

  a VARMABootstrapResult object

- from_index:

  source variable (1-based)

- to_index:

  destination variable (1-based)

- from_lag:

  lag of the source (default 0); must not exceed the AR order p

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

Only the instantaneous and AR (psi) blocks enter the time-expanded
graph; the MA (omega) blocks describe effects of past unobserved
disturbances, which are not nodes of the variable graph (the Python
reference does the same).

## Examples

``` r
s <- generate_varmalingam_sample(n = 300, seed = 42)
bs <- lingam_varma_bootstrap(s$data,
  n_sampling = 5L, order = c(1, 1), criterion = NULL,
  reg_method = "ols", prune = FALSE, seed = 1, verbose = FALSE
)
get_varma_paths(bs, from_index = 1, to_index = 3)
#>      path      effect probability
#> 1 1, 2, 3 -0.26152591           1
#> 2    1, 3 -0.09018331           1
```
