# Select the VARMA order by information criterion

All candidate orders `(p, q)` in `0:max_p x 0:max_q` (excluding
`(0, 0)`) are compared on a **common sample**: the stage-1 long VAR is
estimated once with `h` derived from `(max_p, max_q)`, and every stage-2
candidate regression uses the same window t = h + max_q + 1 .. n, so the
criteria are comparable across orders (same rationale as
[`select_var_lag()`](https://morimotoosamu.github.io/lingamr/reference/select_var_lag.md)).

## Usage

``` r
select_varma_order(X, max_p, max_q, criterion = "bic")
```

## Arguments

- X:

  numeric matrix (n_samples x n_features)

- max_p:

  maximum AR order to consider

- max_q:

  maximum MA order to consider

- criterion:

  "bic", "aic", or "hqic"

## Value

the selected order `c(p, q)` (integer vector)
