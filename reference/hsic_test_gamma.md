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

  numeric vector

- Y:

  numeric vector (same length as X)

## Value

list(stat = HSIC test statistic, p = gamma-approximated p-value)

## Details

The gamma-approximation variance estimator is only defined for n \>= 6
(its closed form divides by `(n-1)(n-2)(n-3)`); below that this errors
instead of silently returning a NaN p-value that would otherwise
propagate into `NA`-valued rejection decisions in callers. A constant
input (zero variance) carries no dependence information, so it is
treated as trivially independent (`p = 1`) rather than routed through
the degenerate kernel-width / centered-Gram-matrix computation.
