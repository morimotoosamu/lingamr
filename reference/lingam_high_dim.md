# High-Dimensional Direct LiNGAM

A variant of Direct LiNGAM for high-dimensional data (large `p`, or
`p > n`). Causal order search is based on moment statistics of
non-Gaussianity rather than pairwise independence measures, which is
considerably faster for many variables. Unlike
[`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md),
the algorithm is deterministic (no random restarts).

## Usage

``` r
lingam_high_dim(X, J = 3L, K = 4L, alpha = 0.5, estimate_adj_mat = TRUE)
```

## Arguments

- X:

  Numeric matrix (n_samples x n_features), data frame or matrix

- J:

  Assumed largest in-degree (single integer, must be \>= 3)

- K:

  Degree of the moment used to measure non-Gaussianity (single integer,
  must be \>= 1)

- alpha:

  Cutoff coefficient for pruning away false parents (single numeric
  value in `[0, 1]`)

- estimate_adj_mat:

  Whether to estimate the adjacency matrix (single logical value). If
  `FALSE`, causal-order search still runs but `adjacency_matrix` is
  returned as an NA-filled matrix (not `NULL`, so downstream S3 methods
  keep working).

## Value

A `LingamResult` object (list), the same class returned by
[`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md),
containing:

- `adjacency_matrix`: adjacency matrix B (n_features x n_features).
  **Convention: `B[i, j]` is the causal coefficient from variable j to
  variable i (j -\> i).** Zero elements indicate no causal relationship.
  All-`NA` when `estimate_adj_mat = FALSE`.

- `causal_order`: estimated causal order (integer vector of 1-based
  indices). Earlier elements are more upstream.

## Details

When `n_samples <= n_features`, the adjacency matrix cannot be estimated
with the usual BIC-based Adaptive LASSO, so this function falls back to
a cross-validated LASSO
([`glmnet::cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html))
after emitting a warning. The upstream Python implementation uses
`LassoLarsCV` for this fallback; `cv.glmnet` follows the same
cross-validation design but is not numerically identical (different
solver: coordinate descent vs. LARS).

The pruning statistic used during causal-order search (Wang & Drton
2020) is the minimum of a conditional non-Gaussianity statistic over
every size-appropriate conditioning subset. The upstream Python
implementation (cdt15/lingam) has a `return` statement mis-indented
inside the loop over conditioning subsets, causing it to only ever
evaluate the first subset. This R implementation intentionally does not
replicate that bug and evaluates every subset, so `causal_order` and
`adjacency_matrix` are not numerically identical to the upstream Python
package.

## References

Wang, Y. S. and Drton, M. (2020). High-dimensional causal discovery
under non-Gaussianity. Biometrika, 107(1), 41-59.

## Examples

``` r
sample <- generate_lingam_sample_6(n = 300, seed = 1)
result <- lingam_high_dim(sample$data)
result$causal_order
#> [1] 4 3 1 5 2 6
round(result$adjacency_matrix, 3)
#>       x0 x1     x2    x3 x4 x5
#> x0 0.000  0  0.000 2.918  0  0
#> x1 2.970  0  2.015 0.000  0  0
#> x2 0.000  0  0.000 5.991  0  0
#> x3 0.000  0  0.000 0.000  0  0
#> x4 8.037  0 -1.001 0.000  0  0
#> x5 4.021  0  0.000 0.000  0  0

# \donttest{
if (requireNamespace("glmnet", quietly = TRUE)) {
  # n <= p: falls back to cross-validated LASSO with a warning
  wide_sample <- generate_lingam_large_sample(p = 30, n = 25, seed = 1)
  wide_result <- lingam_high_dim(wide_sample$data)
  wide_result$causal_order
}
#> Warning: Since n_samples <= n_features, the adjacency matrix is estimated with cross-validated lasso (cv.glmnet) instead of BIC-based lambda selection.
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#>  [1]  2  1  3  8 20 13 27  4  5 12 26  6  7 16  9 15 14 23 10 18 21 11 22 24 29
#> [26] 30 25 19 17 28
# }
```
