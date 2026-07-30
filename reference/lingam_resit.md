# RESIT causal discovery for nonlinear additive noise models

R port of the RESIT algorithm (Regression with Subsequent Independence
Test) from the Python `lingam` package. RESIT assumes a nonlinear
additive noise model `x_i = f_i(parents(x_i)) + e_i` and recovers the
causal structure in two phases: (1) the causal order is estimated by
repeatedly detaching the most sink-like variable, i.e. the variable
whose regression residual on all remaining variables is least dependent
on them (smallest HSIC statistic); (2) superfluous parents are pruned by
testing, for each candidate parent, whether the residual regressed on
the other parents is already independent of the parent set (HSIC p-value
greater than `alpha`).

## Usage

``` r
lingam_resit(X, regressor = "gam", alpha = 0.01, prior_knowledge = NULL)
```

## Arguments

- X:

  numeric matrix or data frame of observed variables

- regressor:

  nonlinear regressor used for all internal regressions. Either the
  string `"gam"` (default; requires the suggested package mgcv, and fits
  a smoothing-spline GAM per regression) or a function `function(X, y)`
  that receives a predictor matrix and a response vector and returns the
  fitted values as a numeric vector of length `nrow(X)` (the model is
  only ever evaluated on its own training data, mirroring
  `regressor.fit(X, y); regressor.predict(X)` in Python)

- alpha:

  significance level of the HSIC independence test used for edge pruning
  (default: 0.01; must be non-negative)

- prior_knowledge:

  optional prior-knowledge matrix with elements 1 (directed path
  exists), 0 (no directed path), and -1 or NA (unknown), with the same
  `[to, from]` orientation as the adjacency matrix

## Value

An object of class `ResitResult` with elements:

- `adjacency_matrix`: (p x p) 0/1 matrix; `B[i, j] = 1` means an edge
  `j -> i` (row = to, col = from). Entries are edge indicators, not
  coefficients.

- `causal_order`: estimated causal order (1-based column positions,
  source first).

- `regressor`: label of the regressor used (`"gam"` or
  `"user function"`).

## Details

Because the model is nonlinear, the returned `adjacency_matrix` contains
0/1 edge indicators, **not** connection strengths, and total causal
effects are undefined: the Python implementation's
[`estimate_total_effect()`](https://morimotoosamu.github.io/lingamr/reference/estimate_total_effect.md)
(always 0) and
[`get_error_independence_p_values()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values.md)
(always a zero matrix) are intentionally not ported.

Note the direction of `alpha`: a parent is removed when the HSIC p-value
exceeds `alpha`, so *larger* values of `alpha` make the test stricter
about declaring independence and therefore *keep more edges*.

The HSIC test in phase 1 measures dependence between a residual and the
joint (multivariate) kernel of up to `ncol(X) - 1` predictors. Following
the Python original, variables are not standardized beforehand; if the
variable scales differ wildly, the largest-scale variable dominates the
kernel distances. Each HSIC call builds n x n Gram matrices, and
`O(ncol(X)^2)` regressions and HSIC tests are performed overall, so the
method is not recommended for `nrow(X)` in the thousands.

At least 6 observations are required (the lower bound of the HSIC
gamma-approximation test), which is stricter than the linear methods in
this package.

Deviation from the Python original when `prior_knowledge` excludes every
remaining sink candidate: the candidate set falls back to all remaining
variables (the original would fail), matching
[`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md).

## References

J. Peters, J. M. Mooij, D. Janzing, B. Schoelkopf. Causal discovery with
continuous additive noise models. Journal of Machine Learning Research,
15: 2009-2053, 2014.

## Examples

``` r
# \donttest{
if (requireNamespace("mgcv", quietly = TRUE)) {
  nonlinear <- generate_resit_sample(n = 300, seed = 1)
  result <- lingam_resit(nonlinear$data)
  print(result)
}
#> RESIT Result
#>   Variables : 4
#>   Regressor : gam
#>   Causal order: x0 -> x1 -> x2 -> x3
#> 
#> Adjacency matrix (row = to, col = from):
#>   (entries are 0/1 edge indicators, not coefficients)
#>    x0 x1 x2 x3
#> x0  0  0  0  0
#> x1  1  0  0  0
#> x2  1  1  0  0
#> x3  0  0  1  0
# }
```
