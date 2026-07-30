# Built-in GAM regressor for RESIT

Fits `mgcv::gam(y ~ s(V1) + s(V2) + ...)` and returns the fitted values.
The basis dimension `k` is capped so that the model stays estimable on
small or resampled data: bootstrap resamples duplicate rows, which
lowers the number of unique covariate values, and mgcv errors out when a
basis has more dimensions than unique values (or more total coefficients
than observations). When even `k = 3` is not affordable the smooth terms
degrade to plain linear terms.

## Usage

``` r
resit_gam_fitted(X, y)
```

## Arguments

- X:

  numeric matrix of predictors (n x d, d \>= 1)

- y:

  numeric response vector (length n)

## Value

fitted values (length n)
