# Special Data: Mixed Variables, Multiple Groups, and Missing Values

Direct LiNGAM assumes a single, complete data matrix of continuous
variables. This article covers three extensions for data that break
those assumptions:

- **LiM**
  ([`lingam_lim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_lim.md)):
  data mixing continuous and discrete variables (binary or Poisson
  counts).
- **MultiGroup Direct LiNGAM**
  ([`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md)):
  several datasets that share a causal structure but not coefficient
  values.
- **Missing data**
  ([`bootstrap_with_imputation()`](https://morimotoosamu.github.io/lingamr/reference/bootstrap_with_imputation.md)):
  bootstrap combined with multiple imputation for data containing `NA`.

``` r

library(lingamr)
```

## LiNGAM for Mixed Data (LiM)

Direct LiNGAM assumes every variable is continuous.
[`lingam_lim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_lim.md)
relaxes this assumption and estimates a causal structure from data
containing a mixture of continuous and discrete variables, following
Zeng et al. (2022). It combines a NOTEARS-style continuous optimization
(the “global” phase) with a combinatorial local search over edge
directions, pruning, and edge addition (the “local” phase).

[`generate_lim_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_lim_sample.md)
generates a small dataset with a known causal chain of continuous and
discrete variables: `x1` (continuous) -\> `x2` (discrete) -\> `x3`
(continuous).

``` r

set.seed(1)
lim_dat <- generate_lim_sample(n = 2000)
head(lim_dat$data)
#>           x1 x2         x3
#> 1  0.1182559  0 -0.8936636
#> 2 -1.8695490  0 -1.2651618
#> 3 -3.3867259  0  2.1530815
#> 4 -0.4395899  1  0.6618645
#> 5  0.3215812  0 -0.4652433
#> 6  1.6721779  0  0.8429619
lim_dat$is_continuous
#> [1]  TRUE FALSE  TRUE
```

[`lingam_lim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_lim.md)
requires `is_continuous`, a logical vector marking which columns are
continuous (`TRUE`) versus discrete (`FALSE`). Because the optimization
starts from a random initial point, reproducibility requires
[`set.seed()`](https://rdrr.io/r/base/Random.html).

``` r

lim_result <- lingam_lim(lim_dat$data, is_continuous = lim_dat$is_continuous)
print(lim_result)
#> LiM Result
#>   Variables : 3
#>   Variable types: continuous, discrete, continuous
#>   Causal order: x1 -> x2 -> x3
#> 
#> Adjacency matrix (row = to, col = from):
#>    x1    x2 x3
#> x1  0 0.000  0
#> x2  1 0.000  0
#> x3  0 1.657  0
```

As with
[`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md),
`adjacency_matrix` follows the `B[i, j]` = j -\> i convention (row = to,
column = from), and `causal_order` lists the estimated topological order
as 1-based indices.

``` r

colnames(lim_dat$data)[lim_result$causal_order]
#> [1] "x1" "x2" "x3"
```

### Count Variables (Poisson)

Discrete variables are binary (0/1) by default; setting
`is_poisson = TRUE` treats them as Poisson-distributed counts instead
(the local search phase then scores them with Poisson regression
log-likelihoods). `generate_lim_sample(is_poisson = TRUE)` generates a
matching count-data example:

``` r

set.seed(2)
lim_pois <- generate_lim_sample(n = 2000, is_poisson = TRUE)
head(lim_pois$data)
#>            x1 x2         x3
#> 1  0.95924714  2  1.9627734
#> 2 -0.27872569  0 -0.3284517
#> 3 -1.63830626  0  0.4878718
#> 4 -0.20204595  1  1.7284788
#> 5  0.02041924  5  3.2550128
#> 6 -0.42685945  1  0.5478110

lim_pois_result <- lingam_lim(lim_pois$data,
  is_continuous = lim_pois$is_continuous, is_poisson = TRUE
)
print(lim_pois_result)
#> LiM Result
#>   Variables : 3
#>   Variable types: continuous, discrete (count), continuous
#>   Causal order: x1 -> x2 -> x3
#> 
#> Adjacency matrix (row = to, col = from):
#>    x1    x2 x3
#> x1  0 0.000  0
#> x2  1 0.000  0
#> x3  0 0.505  0
```

See
[`?lingam_lim`](https://morimotoosamu.github.io/lingamr/reference/lingam_lim.md)
for details on the local phase’s edge-weight convention and its numeric
differences from the Python implementation.

## Multi-Group Direct LiNGAM

[`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
fits a single dataset. When data comes from several sources that
plausibly share the same causal structure but not the same strength of
effect (e.g. the same study run at multiple sites, or the same process
observed in different time periods),
[`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md)
jointly estimates a **common causal order** across all groups while
still allowing each group its own adjacency matrix (structural
coefficients), following Shimizu (2012).

[`generate_multi_group_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_multi_group_sample.md)
generates two datasets that share the causal structure of
[`generate_lingam_sample_6()`](https://morimotoosamu.github.io/lingamr/reference/generate_lingam_sample_6.md)
but with slightly different coefficients per group.

``` r

mg <- generate_multi_group_sample(n = c(1000, 1000), seed = 42)
lapply(mg$data_list, head, 3)
#> $group1
#>         x0        x1       x2        x3        x4        x5
#> 1 2.814924 18.017120 4.543655 0.6333728 18.160090 12.236660
#> 2 1.889685 10.956005 2.188091 0.3175366 13.172754  7.932657
#> 3 1.008905  6.990652 1.953131 0.2409218  6.702107  4.797122
#> 
#> $group2
#>          x0        x1        x2         x3        x4        x5
#> 1 0.7259014  5.482225 0.9301592 0.01259095  5.275903  3.321061
#> 2 2.4321051 17.252303 3.3989705 0.41696287 16.459975 11.368701
#> 3 1.5550457 10.342355 1.8713591 0.24518297 10.520165  7.859683
```

``` r

mg_result <- lingam_multi_group(mg$data_list, reg_method = "ols")
print(mg_result)
#> Multi-Group Direct LiNGAM Result
#>   Groups      : 2 (group1, group2)
#>   Variables   : 6
#>   Causal order (common): x3 -> x0 -> x5 -> x2 -> x4 -> x1
#> 
#> [group1] Adjacency matrix (row = to, col = from):
#>        x0 x1     x2     x3     x4    x5
#> x0  0.000  0  0.000  3.033  0.000 0.000
#> x1  3.237  0  1.965  0.014 -0.034 0.006
#> x2 -0.236  0  0.000  6.112  0.000 0.049
#> x3  0.000  0  0.000  0.000  0.000 0.000
#> x4  7.921  0 -1.063  0.399  0.000 0.018
#> x5  4.016  0  0.000 -0.003  0.000 0.000
#> 
#> [group2] Adjacency matrix (row = to, col = from):
#>       x0 x1     x2     x3    x4     x5
#> x0 0.000  0  0.000  3.504 0.000  0.000
#> x1 2.732  0  2.568  0.083 0.034  0.093
#> x2 0.154  0  0.000  6.322 0.000 -0.024
#> x3 0.000  0  0.000  0.000 0.000  0.000
#> x4 8.483  0 -1.487 -0.110 0.000  0.006
#> x5 4.515  0  0.000 -0.045 0.000  0.000
```

`causal_order` is shared by all groups; `adjacency_matrices` holds one
matrix per group, each following the usual `B[i, j]` = j -\> i
convention.

To analyze a single group with the rest of `lingamr`’s single-group
tooling (total causal effects, independence tests, plotting), extract it
as a plain `LingamResult` with
[`get_group_result()`](https://morimotoosamu.github.io/lingamr/reference/get_group_result.md):

``` r

g1 <- get_group_result(mg_result, "group1")
class(g1)
#> [1] "LingamResult"

estimate_all_total_effects(mg$data_list$group1, g1, method = "ols")
#>             x0 x1        x2        x3          x4          x5
#> x0  0.00000000  0  0.000000  3.033460  0.00000000  0.00000000
#> x1  2.90911952  0  2.001580 21.058733 -0.03397056  0.10299386
#> x2 -0.03933572  0  0.000000  5.992677  0.00000000  0.04894766
#> x3  0.00000000  0  0.000000  0.000000  0.00000000  0.00000000
#> x4  8.03407606  0 -1.062516 18.276121  0.00000000 -0.03416285
#> x5  4.01586857  0  0.000000 12.179395  0.00000000  0.00000000
```

[`lingam_multi_group_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group_bootstrap.md)
provides bootstrap stability estimates in the same joint fashion: every
iteration resamples each group independently, then jointly re-estimates
the causal order and per-group adjacency matrices. It returns a named
list of per-group `BootstrapResult` objects, so the existing bootstrap
query functions apply directly per group:

``` r

mg_bs <- lingam_multi_group_bootstrap(mg$data_list,
  n_sampling = 20L, reg_method = "ols", seed = 1, verbose = FALSE
)
get_probabilities(mg_bs$group1)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]  0.0  0.0 0.30    1 0.00 0.00
#> [2,]  1.0  0.0 1.00    1 0.70 0.80
#> [3,]  0.7  0.0 0.00    1 0.00 0.45
#> [4,]  0.0  0.0 0.00    0 0.00 0.00
#> [5,]  1.0  0.3 1.00    1 0.00 0.55
#> [6,]  1.0  0.2 0.55    1 0.45 0.00
```

Note that
[`lingam_multi_group_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group_bootstrap.md)’s
total causal effects are computed as path-coefficient products over each
iteration’s adjacency matrix, not via regression; this matches the
upstream Python implementation but differs from
[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md)’s
regression-based
[`estimate_all_total_effects()`](https://morimotoosamu.github.io/lingamr/reference/estimate_all_total_effects.md).

## Causal Discovery with Missing Data

All algorithms above assume a complete data matrix. When `X` contains
missing values (`NA`),
[`bootstrap_with_imputation()`](https://morimotoosamu.github.io/lingamr/reference/bootstrap_with_imputation.md)
combines bootstrap resampling with multiple imputation: each resample is
imputed into several complete datasets, and a common causal structure is
jointly estimated across them with
[`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md)
(treating the imputed copies as “groups” that share one causal order).
This is an R port of the Python
`lingam.tools.bootstrap_with_imputation()`.

``` r

sample6_na <- generate_lingam_sample_6(n = 1000, seed = 1)
X_na <- sample6_na$data
set.seed(1)
X_na$x5[sample.int(nrow(X_na), size = round(0.1 * nrow(X_na)))] <- NA # MCAR 10% on x5
```

``` r

bwi <- bootstrap_with_imputation(X_na,
  n_sampling = 20L, n_repeats = 5L, seed = 42, verbose = FALSE
)
print(bwi)
#> ImputationBootstrapResult: 20 samplings x 5 repeats, 6 features, 100 missing cells (original data)
```

The default imputer is `mice::mice(method = "norm")` (Bayesian linear
regression), the closest standard R equivalent of the upstream Python
default (`IterativeImputer(sample_posterior = TRUE)`); numeric results
will not match the Python implementation. Both the imputer and the
causal-discovery step can be swapped for a custom `function` via the
`imputer` and `cd_fit` arguments.

Because each iteration produces `n_repeats` adjacency matrices (one per
imputed dataset), the result’s shape differs from
[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md).
[`as_bootstrap_result()`](https://morimotoosamu.github.io/lingamr/reference/as_bootstrap_result.md)
collapses the `n_repeats` dimension (median or mean) into a regular
`BootstrapResult`, so the existing bootstrap query functions apply as
usual:

``` r

bs_na <- as_bootstrap_result(bwi, aggregate = "median")
get_probabilities(bs_na)
#>    x0 x1 x2 x3 x4 x5
#> x0  0  0  0  1  0  0
#> x1  1  0  1  0  0  0
#> x2  0  0  0  1  0  0
#> x3  0  0  0  0  0  0
#> x4  1  0  1  0  0  0
#> x5  1  0  0  0  0  0
```

[`get_total_causal_effects()`](https://morimotoosamu.github.io/lingamr/reference/get_total_causal_effects.md)
is not available on this `BootstrapResult`, since
[`bootstrap_with_imputation()`](https://morimotoosamu.github.io/lingamr/reference/bootstrap_with_imputation.md)
never computes total effects.

## Related Articles

- [Method selection
  guide](https://morimotoosamu.github.io/lingamr/articles/method-selection.md)
  — which method fits your data
- [Direct LiNGAM in
  depth](https://morimotoosamu.github.io/lingamr/articles/direct-lingam.md)
- [Bootstrap and
  diagnostics](https://morimotoosamu.github.io/lingamr/articles/bootstrap-diagnostics.md)
