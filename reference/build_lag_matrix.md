# Build the lagged design matrix for VAR models

Constructs Z such that column block k (columns `(k-1)*p+1` to `k*p`)
contains `X_{t-k}` for `t = lags+1, ..., n`.

## Usage

``` r
build_lag_matrix(X, lags)
```

## Arguments

- X:

  numeric matrix (n_samples x n_features)

- lags:

  lag order

## Value

matrix of shape `(n - lags, lags * p)`
