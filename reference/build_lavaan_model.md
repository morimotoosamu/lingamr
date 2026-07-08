# Build a lavaan model string from an adjacency matrix

Converts a lingamr-convention adjacency matrix (`B[i, j]` = causal
coefficient from j to i) into a lavaan model syntax string. Non-zero
elements become regression paths (`xi ~ xj`); `NA` elements (used by
e.g.
[`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)
to mark a latent confounder between two variables) become a residual
covariance (`xi ~~ xj`) between the two machine-named variables, which
is the standard lavaan idiom equivalent to a two-indicator latent common
cause with one loading fixed (as used by the Python `semopy`-based
original).

## Usage

``` r
build_lavaan_model(B, var_names)
```

## Arguments

- B:

  adjacency matrix (machine names `x0, x1, ...` expected as row/column
  indices; the caller supplies names via `var_names`)

- var_names:

  machine variable names, length `ncol(B)`

## Value

character scalar, lavaan model syntax (possibly empty string)
