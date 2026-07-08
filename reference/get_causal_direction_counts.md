# Get counts, proportions, and causal effects of causal directions

Get counts, proportions, and causal effects of causal directions

## Usage

``` r
get_causal_direction_counts(
  result,
  n_directions = NULL,
  min_causal_effect = NULL,
  split_by_causal_effect_sign = FALSE,
  labels = NULL
)
```

## Arguments

- result:

  BootstrapResult object

- n_directions:

  How many of the top entries to return (NULL = all)

- min_causal_effect:

  Minimum threshold for the causal effect (NULL = 0)

- split_by_causal_effect_sign:

  Whether to split by the sign of the causal effect

- labels:

  Vector of variable names (NULL allowed; if provided, adds from_name
  and to_name columns)

## Value

A data frame containing the following columns:

- `from`, `to`: 1-based indices of the causal (from) and effect (to)
  variables.

- `count`: Number of bootstrap samples in which this specific causal
  direction was identified.

- `proportion`: The frequency of the direction's occurrence (count /
  n_sampling), representing its bootstrap probability.

- `mean_effect`: The average value of the estimated causal effects
  across samples where this direction was identified.

- `median_effect`: The median value of the estimated causal effects,
  providing a robust estimate of the effect size.

- `sd_effect`: The standard deviation of the causal effect estimates,
  indicating the stability of the effect size.

- `ci_lower`, `ci_upper`: The lower (2.5%) and upper (97.5%) bounds of
  the bootstrap confidence interval for the causal effect.

- `sign` (optional): The sign of the causal effect (1 for positive, -1
  for negative), included if `split_by_causal_effect_sign = TRUE`.

- `from_name`, `to_name` (optional): Character labels for the variables,
  included if `labels` were provided.

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

bs_model <- lingam_direct_bootstrap(LiNGAM_sample_1000$data,
  n_sampling = 30L, reg_method = "ols", seed = 42
)
#> Bootstrap: 30 iterations, method=ols (sequential)
#>   iteration 1 / 30
#>   iteration 10 / 30
#>   iteration 20 / 30
#>   iteration 30 / 30
#> Completed in 0.2 seconds.

get_causal_direction_counts(bs_model, labels = names(LiNGAM_sample_1000$data))
#>    from to count proportion   mean_effect median_effect  sd_effect     ci_lower
#> 1     1  6    30 1.00000000  3.9169074323   4.003899696 0.21847967  3.445944525
#> 2     1  2    29 0.96666667  3.1965221742   3.068221838 0.30992203  2.841587361
#> 3     1  5    29 0.96666667  8.0389161985   8.006451065 0.10854439  7.880941612
#> 4     3  2    29 0.96666667  1.9673174150   1.966453543 0.05027476  1.879777124
#> 5     3  5    29 0.96666667 -1.0480749405  -1.053357169 0.06179135 -1.137049836
#> 6     4  1    29 0.96666667  3.1441078645   3.091442871 0.17659155  2.932002923
#> 7     4  2    29 0.96666667  0.0185014076  -0.008653801 0.22438674 -0.443300570
#> 8     4  3    29 0.96666667  6.0253099490   6.006990357 0.08010239  5.901857344
#> 9     4  5    29 0.96666667  0.4045350222   0.382649748 0.23513863  0.037881264
#> 10    4  6    29 0.96666667 -0.1622345285  -0.079885472 0.24550716 -0.646327934
#> 11    6  2    19 0.63333333  0.0153784233   0.008629175 0.03166615 -0.024662694
#> 12    5  2    18 0.60000000 -0.0500280758  -0.050117494 0.03396618 -0.097710792
#> 13    3  6    17 0.56666667  0.0857324665   0.085303740 0.07256644 -0.052388557
#> 14    6  5    17 0.56666667 -0.0006059292  -0.001975088 0.02392258 -0.041830447
#> 15    3  1    16 0.53333333 -0.0321416604  -0.035631751 0.03269372 -0.087898264
#> 16    1  3    14 0.46666667 -0.2478974798  -0.262812042 0.11310103 -0.438743808
#> 17    5  6    13 0.43333333  0.0349666609   0.029746764 0.03185837 -0.019204672
#> 18    6  3    13 0.43333333  0.0619274966   0.067423011 0.02340254  0.023069303
#> 19    2  5    12 0.40000000  0.0492728268  -0.047501331 0.27277363 -0.064307429
#> 20    2  6    11 0.36666667 -0.0128677846  -0.016891773 0.03405618 -0.061464005
#> 21    1  4     1 0.03333333  0.0111159044   0.011115904 0.00000000  0.011115904
#> 22    2  1     1 0.03333333  0.0428822972   0.042882297 0.00000000  0.042882297
#> 23    2  3     1 0.03333333  0.4080999577   0.408099958 0.00000000  0.408099958
#> 24    2  4     1 0.03333333  0.0010149466   0.001014947 0.00000000  0.001014947
#> 25    3  4     1 0.03333333  0.1360229125   0.136022912 0.00000000  0.136022912
#> 26    5  1     1 0.03333333  0.1075834029   0.107583403 0.00000000  0.107583403
#> 27    5  3     1 0.03333333 -0.1429104530  -0.142910453 0.00000000 -0.142910453
#> 28    5  4     1 0.03333333  0.0094097550   0.009409755 0.00000000  0.009409755
#> 29    6  4     1 0.03333333 -0.0049963432  -0.004996343 0.00000000 -0.004996343
#>        ci_upper from_name to_name
#> 1   4.220378958        x0      x5
#> 2   3.779718648        x0      x1
#> 3   8.208283705        x0      x4
#> 4   2.047119882        x2      x1
#> 5  -0.932471742        x2      x4
#> 6   3.530778717        x3      x0
#> 7   0.357810619        x3      x1
#> 8   6.208911907        x3      x2
#> 9   0.902335801        x3      x4
#> 10  0.165966658        x3      x5
#> 11  0.079210779        x5      x1
#> 12  0.009318859        x4      x1
#> 13  0.190346233        x2      x5
#> 14  0.037238724        x5      x4
#> 15  0.023337357        x2      x0
#> 16 -0.066736406        x0      x2
#> 17  0.080052261        x4      x5
#> 18  0.094281792        x5      x2
#> 19  0.669829251        x1      x4
#> 20  0.051706273        x1      x5
#> 21  0.011115904        x0      x3
#> 22  0.042882297        x1      x0
#> 23  0.408099958        x1      x2
#> 24  0.001014947        x1      x3
#> 25  0.136022912        x2      x3
#> 26  0.107583403        x4      x0
#> 27 -0.142910453        x4      x2
#> 28  0.009409755        x4      x3
#> 29 -0.004996343        x5      x3
```
