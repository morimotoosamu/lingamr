# Build the joined lagged design for [`varma_total_effect_core()`](https://morimotoosamu.github.io/lingamr/reference/varma_total_effect_core.md)

`[X_t, ..., X_{t-(p+from_lag)}, e_{t-1}, ..., e_{t-(q+from_lag)}]`.

## Usage

``` r
varma_joined_design(X, ee_full, order, from_lag)
```

## Arguments

- X:

  numeric matrix (n x m), rows ordered in time

- ee_full:

  LiNGAM residuals, full length (n x m)

- order:

  VARMA order c(p, q)

- from_lag:

  lag of the source variable (non-negative integer)

## Value

matrix (n x m(1 + p + q + 2\*from_lag))
