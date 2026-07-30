# CAM-UV (Causal Additive Models with Unobserved Variables)

A causal discovery method for **nonlinear additive** models that allows
unobserved variables. CAM-UV assumes each observed variable is a
generalized additive function of its parents plus independent noise, and
the causal structure is a DAG. It identifies each variable's direct
parents by scanning variable subsets and testing whether the GAM
regression residual is independent of the candidate parents
([`camuv_find_parents()`](https://morimotoosamu.github.io/lingamr/reference/camuv_find_parents.md));
variable pairs whose residuals remain dependent without an identified
edge are reported as connected by an **unobserved causal path** (UCP: a
directed path through an unobserved variable) or an **unobserved
backdoor path** (UBP: a common unobserved ancestor), and are marked `NA`
in the adjacency matrix rather than oriented.

## Usage

``` r
lingam_camuv(
  X,
  alpha = 0.01,
  num_explanatory_vals = 2L,
  independence = "hsic",
  ind_corr = 0.5,
  prior_knowledge = NULL,
  regressor = "gam"
)
```

## Arguments

- X:

  Numeric matrix (n_samples x n_features), data frame or matrix

- alpha:

  Significance level of the HSIC independence test (used when
  `independence = "hsic"`). Must be non-negative.

- num_explanatory_vals:

  Maximum size of the variable subsets scanned when searching for
  parents (the maximum number of explanatory variables considered
  jointly is `num_explanatory_vals - 1`). Larger values increase
  statistical power but grow combinatorially in cost. Must be an integer
  of 1 or more.

- independence:

  Independence measure: "hsic" (default) uses the HSIC
  gamma-approximation test; "fcorr" uses the F-correlation (kernel
  canonical correlation) and rejects based on `ind_corr` instead of a
  p-value.

- ind_corr:

  Threshold on the F-correlation value, used only when
  `independence = "fcorr"` (independence is declared below this value).
  Must be non-negative. Ignored when `independence = "hsic"`.

- prior_knowledge:

  Optional prior knowledge as variable pairs: a 2-column matrix or a
  list of length-2 vectors, where each pair `c(i, j)` means "variable i
  cannot be a cause of variable j". Indices are **1-based** column
  positions (the Python implementation uses 0-based pairs).

- regressor:

  Nonlinear regressor used for all internal regressions, same interface
  as
  [`lingam_resit()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit.md):
  either the string `"gam"` (default; requires the suggested package
  mgcv) or a function `function(X, y)` returning the fitted values as a
  numeric vector of length `nrow(X)`. The Python implementation
  hardcodes pygam's `LinearGAM`, whose spline basis differs from mgcv's,
  so numerical agreement with Python is not expected (structural results
  are).

## Value

An object of class `CAMUVResult` with elements:

- `adjacency_matrix`: (p x p) matrix. **Convention: `B[i, j] = 1` means
  an edge j -\> i (row = to, col = from)**, same as
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md).
  Entries are 0/1 edge indicators, not coefficients (the causal
  functions are nonlinear). Both entries of a variable pair suspected to
  be connected by a UCP or UBP are `NA`.

- `parents_list`: a list of length `n_features`; element `i` is the
  sorted integer vector of the identified direct parents of variable `i`
  (possibly empty). Like
  [`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md),
  there is no `causal_order`.

- `confounded_pairs`: 2-column integer matrix of the variable pairs
  (1-based column positions) suspected to be connected by a UCP or UBP;
  0 rows if none.

- `regressor`: label of the regressor used (`"gam"` or
  `"user function"`).

## Details

The three stages mirror the Python implementation: (1)
[`camuv_get_neighborhoods()`](https://morimotoosamu.github.io/lingamr/reference/camuv_get_neighborhoods.md)
records which raw variable pairs are dependent; (2)
[`camuv_find_parents()`](https://morimotoosamu.github.io/lingamr/reference/camuv_find_parents.md)
scans variable subsets of size 2 up to `num_explanatory_vals`,
identifying the most sink-like variable of each subset and re-testing
until no new parent is found; (3) remaining dependent pairs without an
identified edge are flagged as UCP/UBP pairs.

Every subset test involves GAM regressions and an HSIC test on `n x n`
Gram matrices, so the method is not recommended for `nrow(X)` in the
thousands, and cost grows combinatorially with `num_explanatory_vals`.

`independence = "fcorr"` is only supported with
`num_explanatory_vals = 2` (the default): larger subsets require testing
a residual against several variables jointly, and the F-correlation is
defined for univariate pairs only (the Python implementation breaks on
that combination as well). The HSIC path supports any
`num_explanatory_vals` via the multivariate kernel.

At least 6 observations are required (the lower bound of the HSIC
gamma-approximation test).

Unlike most other estimators in this package there is no bootstrap
variant, matching the Python implementation, and no total-effect
estimation (effects are nonlinear).

## References

Maeda, T. N. and Shimizu, S. (2021). Causal additive models with
unobserved variables. In Proc. Thirty-Seventh Conference on Uncertainty
in Artificial Intelligence (UAI), PMLR 161: 97-106.

## Examples

``` r
# \donttest{
if (requireNamespace("mgcv", quietly = TRUE)) {
  confounded <- generate_camuv_sample(n = 200, seed = 1)
  result <- lingam_camuv(confounded$data)
  print(result)

  # Pairs connected through unobserved variables are left NA
  result$confounded_pairs
}
#> CAM-UV Result
#>   Variables : 6
#>   Regressor : gam
#> 
#> Parent sets:
#>   P(x0) = {}
#>   P(x1) = {x0}
#>   P(x2) = {}
#>   P(x3) = {x0}
#>   P(x4) = {x2}
#>   P(x5) = {}
#> 
#> Pairs with an unobserved causal/backdoor path (UCP/UBP):
#>   x2 -- x5
#>   x3 -- x4
#> 
#> Adjacency matrix (row = to, col = from):
#>   (entries are 0/1 edge indicators, not coefficients;
#>    NA = pair connected through an unobserved variable)
#>    x0 x1 x2 x3 x4 x5
#> x0  0  0  0  0  0  0
#> x1  1  0  0  0  0  0
#> x2  0  0  0  0  0 NA
#> x3  1  0  0  0 NA  0
#> x4  0  0  1 NA  0  0
#> x5  0  0 NA  0  0  0
#>      var1 var2
#> [1,]    3    6
#> [2,]    4    5
# }
```
