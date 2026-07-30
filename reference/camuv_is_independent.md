# Independence judgment at the configured threshold

Faithful port of `CAMUV._is_independent()`: alpha for hsic, ind_corr for
fcorr.

## Usage

``` r
camuv_is_independent(X, Y, independence, alpha, ind_corr)
```

## Arguments

- X:

  numeric vector or single/multi-column matrix

- Y:

  numeric vector or single/multi-column matrix

- independence:

  "hsic" or "fcorr"

- alpha:

  significance level (hsic only)

- ind_corr:

  rejection threshold (fcorr only)

## Value

TRUE if independent
