# RESIT step 2: remove superfluous edges

Faithful port of `RESIT._remove_edges()` (resit.py). For each variable
in causal order (skipping the source) and each of its candidate parents
`l`, the variable is regressed on its current parents excluding `l`; if
the residual is independent (HSIC p-value \> alpha) of the variable's
original parent set, the edge from `l` is dropped.

## Usage

``` r
resit_remove_edges(X, pa, pi_order, reg_fn, alpha)
```

## Arguments

- X:

  data matrix

- pa:

  list of parent index vectors from
  [`resit_estimate_order()`](https://morimotoosamu.github.io/lingamr/reference/resit_estimate_order.md)

- pi_order:

  causal order (source first)

- reg_fn:

  regressor function

- alpha:

  significance level of the HSIC test

## Value

pruned `pa` list

## Details

Note the two distinct parent sets, mirroring the original exactly:
`parents_snapshot` is the parent set frozen before the inner loop and is
always the second HSIC argument, while `predictors` is derived from the
*current* (already pruned) parent set minus `l`.
