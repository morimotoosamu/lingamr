# Compute residuals (error terms) of a LiNGAM model

After validating the inputs (that `lingam_result` is a `LingamResult`,
that X is numeric, and that the dimensions match), returns
`E = X - X B^T`. Shared by the residual-based diagnostic functions.

## Usage

``` r
lingam_residuals(X, lingam_result)
```

## Arguments

- X:

  original data (matrix or data.frame)

- lingam_result:

  return value of
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)

## Value

residual matrix (n_samples x n_features). Preserves the column names of
X.
