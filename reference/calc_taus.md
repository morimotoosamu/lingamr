# Conditional non-Gaussianity statistics (tau), minimized over all conditioning subsets

Per Wang & Drton (2020), the pruning statistic for a candidate variable
is the minimum of the conditional tau statistic over every
size-appropriate conditioning subset. Upstream HighDimDirectLiNGAM
(cdt15/lingam) has a `return` mis-indented inside its loop over
conditioning sets, so it only ever evaluates the first subset; this R
port intentionally does NOT replicate that bug and instead evaluates
every subset in `cond_sets`, so results differ numerically from the
Python package (see `dev/high-dim-direct-lingam-implementation.md`).

## Usage

``` r
calc_taus(Y, yty, pa, ch, k, cond_sets, an_sets)
```

## Arguments

- Y:

  data matrix

- yty:

  cached Gram matrix `t(Y) %*% Y`

- pa:

  index of the candidate parent variable (scalar)

- ch:

  indices of candidate child variables

- k:

  moment degree

- cond_sets:

  list of conditioning sets (each an integer vector of 1-based column
  indices, always including `last_root`)

- an_sets:

  list, parallel to `cond_sets`, of ancestor-candidate sets
  (conditioning-set variables excluded from the corresponding subset)

## Value

numeric vector of length `ncol(Y)`, minimum statistic per variable
