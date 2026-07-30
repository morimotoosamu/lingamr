# Enumerate and aggregate paths over time-expanded bootstrap graphs

Shared core of
[`get_var_paths()`](https://morimotoosamu.github.io/lingamr/reference/get_var_paths.md)
and
[`get_varma_paths()`](https://morimotoosamu.github.io/lingamr/reference/get_varma_paths.md):
builds the time-expanded square graph for each bootstrap adjacency
matrix, enumerates all directed paths from the source (at `from_lag`) to
the destination (at `to_lag`), and aggregates each distinct path's
bootstrap probability and median effect.

## Usage

``` r
collect_time_expanded_paths(
  ams,
  nf,
  n_lags,
  from_index,
  to_index,
  from_lag,
  to_lag,
  min_causal_effect
)
```

## Arguments

- ams:

  list of joined lag matrices (n_features x n_features\*(1 + n_lags));
  block k + 1 holds the lag-k coefficients

- nf:

  number of features

- n_lags:

  number of lag blocks beyond the instantaneous one

- from_index:

  source variable (1-based)

- to_index:

  destination variable (1-based)

- from_lag:

  lag of the source

- to_lag:

  lag of the destination

- min_causal_effect:

  minimum \|effect\| threshold

## Value

a data frame (path, effect, probability), one row per distinct path
