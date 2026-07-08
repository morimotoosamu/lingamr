# Moore-Penrose pseudo-inverse via SVD

Internal replacement for `numpy.linalg.pinv`, used to solve the
conditional regression coefficients from Gram-matrix submatrices in
[`lingam_high_dim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_high_dim.md).
Implemented in base R to avoid a new dependency.

## Usage

``` r
pinv(A, tol = max(dim(A)) * .Machine$double.eps)
```

## Arguments

- A:

  numeric matrix

- tol:

  singular-value cutoff, relative to the largest singular value

## Value

the pseudo-inverse of `A`
