# Select the most sink-like child within a variable set

Faithful port of `CAMUV._get_child()`. For each candidate child, the
child is regressed on the candidate parents plus its already-identified
parents, and the residual is tested against the *residual matrix*
columns of the candidate parents. The candidate whose residual is most
independent wins; `prev_independence` starts at 0 for hsic (p-value
scale, higher is more independent) and at 1 for fcorr (lower is more
independent), and each accepted candidate raises the bar for the next.

## Usage

``` r
camuv_get_child(
  X,
  vars,
  P,
  N,
  Y,
  pk_forbidden,
  independence,
  alpha,
  ind_corr,
  reg_fn,
  get_residual = NULL
)
```

## Arguments

- X:

  data matrix

- vars:

  variable subset (integer vector)

- P:

  current parent list

- N:

  neighborhood list

- Y:

  current residual matrix

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

- get_residual:

  function(v, ids) returning the residual of `X[, v]` regressed on
  `X[, ids]`; defaults to an uncached
  [`camuv_get_residual()`](https://morimotoosamu.github.io/lingamr/reference/camuv_get_residual.md)
  call.
  [`camuv_find_parents()`](https://morimotoosamu.github.io/lingamr/reference/camuv_find_parents.md)
  passes a memoized version, since the same (child, predictor set)
  recurs across subset rescans.

## Value

list(child = index or NA, independent = logical); `independent` reports
whether the winning candidate's dependence value clears the configured
threshold (alpha / ind_corr)
