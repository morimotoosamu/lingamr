# RESIT step 1: estimate the topological (causal) order

Faithful port of `RESIT._estimate_order()` (resit.py). Repeatedly finds
the most sink-like variable: the one whose regression residual on all
other remaining variables has the smallest HSIC dependence statistic.
Sinks are identified first and prepended, so the returned order runs
from source to sink.

## Usage

``` r
resit_estimate_order(X, reg_fn, partial_orders, Aknw)
```

## Arguments

- X:

  data matrix

- reg_fn:

  regressor function from
  [`resit_make_regressor()`](https://morimotoosamu.github.io/lingamr/reference/resit_make_regressor.md)

- partial_orders:

  (from, to) matrix from
  [`extract_partial_orders()`](https://morimotoosamu.github.io/lingamr/reference/extract_partial_orders.md),
  or NULL

- Aknw:

  validated prior-knowledge matrix (negative entries already NA), or
  NULL

## Value

list(pa = list of parent index vectors per variable, order = causal
order, source first)

## Details

Deviation from the Python original:
[`parce_search_candidate()`](https://morimotoosamu.github.io/lingamr/reference/parce_search_candidate.md)
falls back to the full remaining set when prior knowledge excludes every
candidate (the Python code would crash on an empty argmin).
