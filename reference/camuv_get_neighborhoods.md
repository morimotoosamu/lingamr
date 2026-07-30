# Pairwise dependence neighborhoods

Faithful port of `CAMUV._get_neighborhoods()`: `N[[i]]` holds the
variables whose raw column is *dependent* on column i.

## Usage

``` r
camuv_get_neighborhoods(X, independence, alpha, ind_corr)
```

## Arguments

- X:

  data matrix

- independence:

  "hsic" or "fcorr"

- alpha:

  significance level (hsic only)

- ind_corr:

  rejection threshold (fcorr only)

## Value

list of length `ncol(X)` of integer vectors
