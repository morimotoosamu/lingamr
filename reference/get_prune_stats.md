# Enumerate conditioning subsets and compute pruning statistics for one candidate variable

Enumerate conditioning subsets and compute pruning statistics for one
candidate variable

## Usage

``` r
get_prune_stats(Y, yty, i, j, K, last_root, condition_set, J)
```

## Arguments

- Y:

  data matrix

- yty:

  cached Gram matrix `t(Y) %*% Y`

- i:

  candidate variable (scalar, 1-based index)

- j:

  current candidate set (`psi`; `i` is removed internally)

- K:

  moment degree

- last_root:

  most recently fixed causal-order variable, or `NULL` on the first
  iteration

- condition_set:

  conditioning-set variables (already includes `last_root` when
  non-empty)

- J:

  assumed largest in-degree

## Value

numeric vector of length `ncol(Y)`
