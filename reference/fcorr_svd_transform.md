# Low-rank SVD transform used inside [`f_correlation()`](https://morimotoosamu.github.io/lingamr/reference/f_correlation.md)

Given the (already column-centered) incomplete-Cholesky factor `G`,
returns the orthonormalized basis `U` and the shrinkage vector `R` used
to assemble the block canonical-correlation matrix `R_kappa`.

## Usage

``` r
fcorr_svd_transform(G, kappa, n)
```

## Arguments

- G:

  n x d centered incomplete-Cholesky factor

- kappa:

  regularization parameter

- n:

  sample size

## Value

list(U = n x d' matrix, R = length d' vector); `R` has length 0 if no
eigenvalue of `crossprod(G)` clears the `kappa` threshold (i.e. `G`
carries no retainable rank), which callers must handle explicitly.
