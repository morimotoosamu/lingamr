# HSIC gamma test from precomputed parts

Bit-identical to `hsic_test_gamma(X, Y)` when `pre_x`/`pre_y` are
[`hsic_precompute()`](https://morimotoosamu.github.io/lingamr/reference/hsic_precompute.md)
of the same arguments in the same roles. The roles matter at the
last-ulp level: the `mean_` expression is not symmetric in floating
point, so `pre_x` must correspond to the first argument of the original
call being replaced.

## Usage

``` r
hsic_gamma_from_pre(pre_x, pre_y)
```

## Arguments

- pre_x:

  precomputed first argument (from
  [`hsic_precompute()`](https://morimotoosamu.github.io/lingamr/reference/hsic_precompute.md))

- pre_y:

  precomputed second argument

## Value

list(stat = HSIC test statistic, p = gamma-approximated p-value)
