# Reject constant or perfectly collinear columns

Constant and linearly dependent columns produce divisions by zero inside
the pairwise regressions / standardization of the causal order search,
which would otherwise surface as cryptic errors (e.g. "argument is of
length zero") deep in the algorithm. The rank check runs on the centered
matrix so that a column equal to another column plus a constant offset
is also caught.

## Usage

``` r
validate_no_degenerate_columns(X)
```

## Arguments

- X:

  numeric matrix

## Value

`NULL`, invisibly. Stops with an informative error on violation.
