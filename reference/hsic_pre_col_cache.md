# Lazy per-column cache of [`hsic_precompute()`](https://morimotoosamu.github.io/lingamr/reference/hsic_precompute.md) results

Returns a closure `f(k)` that computes `hsic_precompute(X[, k])` on
first use and returns the cached object afterwards. For callers that
test many pairs of columns of a fixed matrix.

## Usage

``` r
hsic_pre_col_cache(X)
```

## Arguments

- X:

  numeric matrix whose columns will be tested

## Value

function(k) returning the precompute object for column k
