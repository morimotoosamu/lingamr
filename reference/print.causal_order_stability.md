# print method for causal_order_stability

print method for causal_order_stability

## Usage

``` r
# S3 method for class 'causal_order_stability'
print(x, ...)
```

## Arguments

- x:

  A `causal_order_stability` object

- ...:

  Additional arguments (unused)

## Value

The input object `x`, invisibly.

## Examples

``` r
dat <- generate_lingam_sample_6()
bs <- lingam_direct_bootstrap(dat$data, n_sampling = 30L, reg_method = "ols", seed = 42)
#> Bootstrap: 30 iterations, method=ols (sequential)
#>   iteration 1 / 30
#>   iteration 10 / 30
#>   iteration 20 / 30
#>   iteration 30 / 30
#> Completed in 0.2 seconds.
print(get_causal_order_stability(bs, labels = names(dat$data)))
#> === Causal Order Stability ===
#> Bootstrap samples:       30
#> Overall stability score: 0.680  (0 = random, 1 = fully stable)
#> 
#> Rank summary (sorted by mean rank; 1 = most upstream):
#>  variable mean_rank sd_rank median_rank mode_rank
#>        x3      1.17    0.91         1.0         1
#>        x0      2.57    0.57         3.0         3
#>        x2      2.93    0.98         2.5         2
#>        x5      4.33    1.32         4.0         3
#>        x4      4.87    0.86         5.0         5
#>        x1      5.13    1.11         5.0         6
#> 
#> Precedence probability P[i, j] = P(variable i precedes j):
#>      x0   x1   x2   x3   x4   x5
#> x0 0.00 0.97 0.47 0.03 0.97 1.00
#> x1 0.03 0.00 0.03 0.03 0.40 0.37
#> x2 0.53 0.97 0.00 0.03 0.97 0.57
#> x3 0.97 0.97 0.97 0.00 0.97 0.97
#> x4 0.03 0.60 0.03 0.03 0.00 0.43
#> x5 0.00 0.63 0.43 0.03 0.57 0.00
```
