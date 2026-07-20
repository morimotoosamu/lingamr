# Convert an ImputationBootstrapResult to a tidy data.frame

Collapses the imputation dimension with
[`as_bootstrap_result()`](https://morimotoosamu.github.io/lingamr/reference/as_bootstrap_result.md)
and then summarizes the causal direction counts like
[`tidy.BootstrapResult()`](https://morimotoosamu.github.io/lingamr/reference/tidy.BootstrapResult.md).

## Usage

``` r
# S3 method for class 'ImputationBootstrapResult'
tidy(x, aggregate = c("median", "mean"), ...)
```

## Arguments

- x:

  The return value of
  [`bootstrap_with_imputation()`](https://morimotoosamu.github.io/lingamr/reference/bootstrap_with_imputation.md)
  (an `ImputationBootstrapResult` object)

- aggregate:

  How to collapse the `n_repeats` imputation dimension, passed to
  [`as_bootstrap_result()`](https://morimotoosamu.github.io/lingamr/reference/as_bootstrap_result.md)
  ("median" or "mean")

- ...:

  Arguments passed to
  [`get_causal_direction_counts()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_direction_counts.md)

## Value

data.frame (from, to, count, proportion, ...)

## Examples

``` r
dat <- generate_lingam_sample_6(n = 200, seed = 1)$data
dat[sample(nrow(dat), 20), 1] <- NA
bs <- bootstrap_with_imputation(dat, n_sampling = 5L, n_repeats = 2L, seed = 42)
#> Bootstrap with imputation: 5 iterations, n_repeats=2 (sequential)
#>   iteration 1 / 5
#> Completed in 0.3 seconds.
tidy(bs)
#>    from to count proportion mean_effect median_effect   sd_effect   ci_lower
#> 1     3  2     5        1.0   2.1383915     2.0617327 0.176654962  1.9992568
#> 2     4  3     5        1.0   5.9970268     5.9873049 0.057209867  5.9427038
#> 3     4  5     5        1.0  18.7973718    18.1917914 1.926593886 17.4957668
#> 4     5  1     5        1.0   0.1552270     0.1565516 0.004672225  0.1479719
#> 5     1  2     4        0.8   2.6377799     2.8992287 0.564333100  1.8736513
#> 6     3  6     4        0.8   0.4542938     0.5195260 0.152648364  0.2481890
#> 7     5  6     4        0.8   0.4893890     0.4897081 0.006038851  0.4832631
#> 8     5  2     2        0.4   0.2463576     0.2463576 0.158645772  0.1397871
#> 9     3  5     1        0.2  -0.7811145    -0.7811145 0.000000000 -0.7811145
#> 10    1  6     1        0.2   3.9956585     3.9956585 0.000000000  3.9956585
#> 11    4  6     1        0.2   1.9828259     1.9828259 0.000000000  1.9828259
#>      ci_upper
#> 1   2.3998647
#> 2   6.0791722
#> 3  21.8022401
#> 4   0.1585205
#> 5   2.9574454
#> 6   0.5495037
#> 7   0.4949723
#> 8   0.3529282
#> 9  -0.7811145
#> 10  3.9956585
#> 11  1.9828259
```
