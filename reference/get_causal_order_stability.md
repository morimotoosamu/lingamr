# Evaluate the stability of the causal order from bootstrap

Aggregates the causal order (causal_order) estimated in each bootstrap
sample and quantifies how stable the order is. Returns the rank
distribution of each variable, the precedence probabilities for variable
pairs, and an overall stability score.

## Usage

``` r
get_causal_order_stability(result, labels = NULL)
```

## Arguments

- result:

  A BootstrapResult object (run with the current version)

- labels:

  A vector of variable names (if NULL, x0, x1, ... are generated
  automatically)

## Value

A list of class `causal_order_stability`, containing:

- `rank_summary`: A summary of the rank of each variable (variable,
  mean_rank, sd_rank, median_rank, mode_rank). Sorted in ascending order
  of mean_rank (from upstream). A rank of 1 is the most upstream.

- `precedence_matrix`: A precedence probability matrix. `P[i, j]` is the
  proportion of bootstrap samples in which variable i was located
  upstream of (before) variable j.

- `stability_score`: An overall stability score, from 0 (random order)
  to 1 (order agrees across all samples). The closer the precedence
  probability of each variable pair is to 0/1, the higher the score.

- `n_sampling`: The number of bootstrap samples.

## Examples

``` r
dat <- generate_lingam_sample_6()
bs <- lingam_direct_bootstrap(dat$data, n_sampling = 30L, seed = 42)
#> Bootstrap: 30 iterations, method=adaptive_lasso (sequential)
#>   iteration 1 / 30
#>   iteration 10 / 30
#>   iteration 20 / 30
#>   iteration 30 / 30
#> Completed in 0.7 seconds.
get_causal_order_stability(bs, labels = names(dat$data))
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
