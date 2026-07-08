# Bottom-up causal order search

Faithful port of `_search_causal_order()` (original 188-225 lines).
Repeatedly finds the most sink-like remaining variable and appends it to
the front of `K_bttm` (bottom-up, so more recently placed variables are
more upstream). Stops as soon as a candidate is rejected by the
independence test; the remaining undetermined variables are returned as
`U_res`.

## Usage

``` r
parce_search_causal_order(
  X,
  U,
  partial_orders,
  independence,
  thresh_p,
  ind_corr
)
```

## Arguments

- X:

  (centered) data matrix

- U:

  all variable indices

- partial_orders:

  matrix of (from, to) pairs, or NULL

- independence:

  "hsic" or "fcorr"

- thresh_p:

  Bonferroni-corrected significance threshold (hsic only)

- ind_corr:

  F-correlation rejection threshold (fcorr only)

## Value

list(K_bttm = integer vector, p_bttm = numeric vector, U_res = integer
vector)
