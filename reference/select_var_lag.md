# Select the VAR lag order by information criterion

All candidate lag orders are compared on a **common sample**: the first
`max_lag` observations are dropped for every candidate so that each
VAR(lag) is estimated over the same time window (t = max_lag + 1 .. n).
This mirrors statsmodels' `VAR.select_order` and makes the criteria
comparable across lags (otherwise a longer lag would be scored on fewer
observations).

## Usage

``` r
select_var_lag(X, max_lag, criterion = "bic")
```

## Arguments

- X:

  numeric matrix (n_samples x n_features)

- max_lag:

  maximum lag order to consider

- criterion:

  "bic", "aic", "hqic", or "fpe"

## Value

the selected lag order (integer)
