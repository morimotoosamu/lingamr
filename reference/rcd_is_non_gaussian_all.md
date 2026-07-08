# Non-Gaussianity judgment (Shapiro-Wilk) for a set of columns

When `n > SHAPIRO_MAX_N`,
[`stats::shapiro.test()`](https://rdrr.io/r/stats/shapiro.test.html)
cannot be used directly (hard cap at 5000), so the deterministic
evenly-spaced subsample from
[`shapiro_subsample()`](https://morimotoosamu.github.io/lingamr/reference/shapiro_subsample.md)
(defined in `R/get_error_independence_p_values.r`) is tested instead,
matching the behavior of
[`test_residual_normality()`](https://morimotoosamu.github.io/lingamr/reference/test_residual_normality.md).
The deterministic thinning keeps results reproducible across calls and
leaves the caller's RNG stream untouched.

## Usage

``` r
rcd_is_non_gaussian_all(Y, cols, shapiro_alpha)
```

## Arguments

- Y:

  data matrix

- cols:

  column indices to test

- shapiro_alpha:

  significance level

## Value

TRUE if all columns reject normality (p \<= shapiro_alpha)
