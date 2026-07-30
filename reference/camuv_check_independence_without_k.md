# Check that the child's residual is dependent on each parent's residual

Faithful port of `CAMUV._check_independence_withou_K()` (typo in the
upstream method name). If the child's current residual is already
independent of some candidate parent's residual *without* regressing on
the candidate set K, that parent adds no information and the candidate
set is rejected.

## Usage

``` r
camuv_check_independence_without_k(
  parents,
  child,
  Y,
  independence,
  alpha,
  ind_corr
)
```

## Arguments

- parents:

  candidate parent indices

- child:

  candidate child index

- Y:

  current residual matrix

- independence:

  "hsic" or "fcorr"

- alpha:

  significance level (hsic only)

- ind_corr:

  rejection threshold (fcorr only)

## Value

TRUE if the child's residual is dependent on every parent's residual
