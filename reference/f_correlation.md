# F-correlation (kernel canonical correlation) between two variables

Bach & Jordan (2002) kernel canonical correlation, as used by
BottomUpParceLiNGAM's `independence = "fcorr"` option. Returns a value
in (roughly) `[0, 1]`; larger means more dependent.

## Usage

``` r
f_correlation(x, y)
```

## Arguments

- x:

  numeric vector

- y:

  numeric vector (same length as x)

## Value

F-correlation value (scalar)

## Details

A constant `x` or `y` (zero variance) carries no dependence information
and is treated as trivially independent (returns `0`) rather than
propagating a division-by-zero standardization into NaN/Inf.
