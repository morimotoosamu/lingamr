# Call a regressor function and validate its return value

User-supplied regressors are called many times deep inside the order
search, so a malformed return value is caught here with a clear message
instead of surfacing as a cryptic arithmetic error.

## Usage

``` r
resit_fitted(reg_fn, Xp, y)
```

## Arguments

- reg_fn:

  function `(X, y) -> fitted values`

- Xp:

  predictor matrix (n x d)

- y:

  response vector (length n)

## Value

fitted values as a plain numeric vector
