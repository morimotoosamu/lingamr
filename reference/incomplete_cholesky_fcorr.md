# Pivoted incomplete Cholesky decomposition of a Gaussian kernel matrix, as used by [`f_correlation()`](https://morimotoosamu.github.io/lingamr/reference/f_correlation.md)

Same greedy pivoted algorithm as
[`incomplete_cholesky_gauss()`](https://morimotoosamu.github.io/lingamr/reference/incomplete_cholesky_gauss.md)
(search_causal_order.r), but with the stopping rule used by the upstream
`_f_correlation.py`: continue while the *sum* of the remaining diagonal
residuals exceeds `tol` (no rank cap), rather than stopping once the
single largest residual drops below a fixed tolerance.

## Usage

``` r
incomplete_cholesky_fcorr(x, sigma, tol)
```

## Arguments

- x:

  input vector (length n)

- sigma:

  width of the Gaussian kernel

- tol:

  stop once the sum of the remaining diagonal residuals is at or below
  this

## Value

n x d matrix `G` with `tcrossprod(G) ~= K`
