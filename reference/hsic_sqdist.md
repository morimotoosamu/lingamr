# Pairwise squared Euclidean distances between rows

The single-column case keeps the original
[`outer()`](https://rdrr.io/r/base/outer.html) formulation so that
existing univariate callers (Parce / RCD) get bit-identical results; the
multivariate case uses the rowSums/tcrossprod expansion, which is only
reached by matrix inputs (new with the RESIT port).

## Usage

``` r
hsic_sqdist(X)
```

## Arguments

- X:

  numeric matrix (n x d)

## Value

n x n matrix of squared distances
