# Independence judgment with an explicit threshold, returning the value

Faithful port of `CAMUV._is_independent_by()`. For `"hsic"` the value is
the gamma-approximated p-value and independence means
`value > threshold`; for `"fcorr"` the value is the F-correlation and
independence means `value < threshold`. Note the opposite directions.

## Usage

``` r
camuv_is_independent_by(X, Y, threshold, independence)
```

## Arguments

- X:

  numeric vector or single/multi-column matrix

- Y:

  numeric vector or single/multi-column matrix

- threshold:

  rejection threshold (see above)

- independence:

  "hsic" or "fcorr"

## Value

list(independent = logical, value = numeric)

## Details

[`f_correlation()`](https://morimotoosamu.github.io/lingamr/reference/f_correlation.md)
is univariate, so the fcorr path requires single-column inputs;
[`lingam_camuv()`](https://morimotoosamu.github.io/lingamr/reference/lingam_camuv.md)
rejects `independence = "fcorr"` with `num_explanatory_vals > 2` up
front, which is the only way a multivariate input could reach this point
(the upstream implementation breaks on that combination too).
