# CAM-UV main search: identify each variable's parents

Faithful port of `CAMUV._find_parents()`. Scans variable subsets of
increasing size `t` (starting at 2); whenever any parent is identified,
`t` resets to 2, otherwise `t` grows until it exceeds
`num_explanatory_vals`. The residual matrix `Y` is updated in place each
time a variable gains a parent. A final pruning pass removes parents
whose residuals turn out to be independent of the child's residual.

## Usage

``` r
camuv_find_parents(
  X,
  maxnum_vals,
  N,
  pk_forbidden,
  independence,
  alpha,
  ind_corr,
  reg_fn
)
```

## Arguments

- X:

  data matrix

- maxnum_vals:

  maximum subset size (`num_explanatory_vals`)

- N:

  neighborhood list from
  [`camuv_get_neighborhoods()`](https://morimotoosamu.github.io/lingamr/reference/camuv_get_neighborhoods.md)

- pk_forbidden:

  forbidden-cause list, or NULL

- independence:

  "hsic" or "fcorr"

- alpha:

  significance level (hsic only)

- ind_corr:

  rejection threshold (fcorr only)

- reg_fn:

  regressor function

## Value

list of parent index vectors per variable
