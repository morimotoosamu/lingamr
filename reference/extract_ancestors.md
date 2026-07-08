# Extract the ancestor sets M(x_i) for every variable

Faithful port of `rcd.py`'s `extract_ancestors()` (original 254-316
lines). Repeatedly scans variable subsets of increasing size, growing
each variable's ancestor set whenever a unique "sink" is found within a
subset. See `dev/rcd-implementation.md` section 3.2 for the full
algorithm description and the reasoning behind the cache-update
position.

## Usage

``` r
extract_ancestors(
  X,
  max_explanatory_num,
  cor_alpha,
  ind_alpha,
  shapiro_alpha,
  MLHSICR,
  independence,
  ind_corr
)
```

## Arguments

- X:

  (uncentered) data matrix

- max_explanatory_num:

  maximum subset size minus 1

- cor_alpha:

  correlation-test significance level

- ind_alpha:

  independence-test significance level (hsic)

- shapiro_alpha:

  non-Gaussianity significance level

- MLHSICR:

  whether to use MLHSICR regression as a fallback

- independence:

  "hsic" or "fcorr"

- ind_corr:

  F-correlation rejection threshold (fcorr)

## Value

list of length `ncol(X)`; each element is a sorted integer vector of
ancestor indices (possibly empty)
