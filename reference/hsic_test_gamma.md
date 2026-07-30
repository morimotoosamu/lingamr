# HSIC independence test with gamma approximation

Faithful port of `hsic_test_gamma()` (hsic.py). O(n^2) in the sample
size because it forms the full n x n Gram matrices; not recommended for
n in the thousands.

## Usage

``` r
hsic_test_gamma(X, Y)
```

## Arguments

- X:

  numeric vector or matrix (n x d)

- Y:

  numeric vector or matrix with the same number of rows as X

## Value

list(stat = HSIC test statistic, p = gamma-approximated p-value)

## Details

Either argument may be a matrix (n x d); its columns are combined into a
single Gaussian kernel over row-wise Euclidean distances, exactly as in
the upstream multivariate implementation (used by RESIT, which tests a
residual against the joint distribution of several predictors). The
variables are not standardized beforehand (upstream behavior), so with
wildly different column scales the largest-scale column dominates the
distance.

The gamma-approximation variance estimator is only defined for n \>= 6
(its closed form divides by `(n-1)(n-2)(n-3)`); below that this errors
instead of silently returning a NaN p-value that would otherwise
propagate into `NA`-valued rejection decisions in callers. A constant
input (zero variance in every column) carries no dependence information,
so it is treated as trivially independent (`p = 1`) rather than routed
through the degenerate kernel-width / centered-Gram-matrix computation.
