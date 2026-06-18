# Build a zero-adjacency LingamResult

A stand-in `LingamResult` whose adjacency matrix is all zeros, so that
[`lingam_residuals()`](https://morimotoosamu.github.io/lingamr/reference/lingam_residuals.md)
returns its input unchanged. This lets the VAR diagnostics reuse the
Direct LiNGAM residual routines on an already-computed residual matrix.

## Usage

``` r
zero_lingam_result(p)
```

## Arguments

- p:

  number of features

## Value

a `LingamResult` with a p x p zero adjacency matrix
