# Latent Confounders: ParceLiNGAM and RCD

Standard Direct LiNGAM assumes there is no latent (unobserved)
confounder: any variable that causes two or more of the observed
variables must itself be observed. When that assumption fails,
[`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
will still return a full causal order, but silently — some part of it
may be wrong, with no indication of which part.

This article covers the two **linear** methods that address latent
confounders:

- **BottomUpParceLiNGAM**
  ([`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)):
  returns the variables it could not order as a single *unresolved
  block*, signaling where a latent confounder is likely.
- **RCD**
  ([`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)):
  estimates each variable’s ancestor set and flags *specific pairs* that
  share a latent confounder.

If the causal relationships may also be **nonlinear**, see CAM-UV in the
[nonlinear methods
article](https://morimotoosamu.github.io/lingamr/articles/nonlinear.md),
which handles unobserved variables in nonlinear models.

``` r

library(lingamr)
```

## Bottom-Up ParceLiNGAM

[`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)
(BottomUpParceLiNGAM, Tashiro et al. 2014) searches for the causal order
from the sink (most downstream) side, testing at each step whether a
candidate variable’s residual is independent of the others. As soon as
that test is rejected, the search stops, and every variable it could not
yet place is returned together as a single **unresolved block** — a
signal that those variables likely share a latent confounder, rather
than a (possibly wrong) guess at their order.

[`generate_parce_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_parce_sample.md)
generates a 7-variable model in which `x6` is an unobserved common cause
of `x2` and `x3`; only `x0`-`x5` are returned as data.

``` r

# HSIC is O(n^2), so a moderate n keeps this article fast to build
confounded <- generate_parce_sample(n = 500, seed = 1)
head(confounded$data)
#>          x0        x1        x2        x3          x4        x5
#> 1 0.6154746 1.7104554 1.0618261 1.0851944  0.57917376 0.8337563
#> 2 1.5905703 2.4770365 1.4291087 1.4325230  0.56022597 0.8686206
#> 3 1.1007549 2.1817888 1.5289901 1.8037643 -0.03671602 1.4001192
#> 4 1.7744689 2.7106515 2.7714036 2.4797583 -0.10133399 1.3102925
#> 5 0.5433612 0.7244786 0.5217204 0.8755981  0.82504738 1.2597767
#> 6 1.7671488 1.8838085 1.8358794 2.7663075  0.90699176 1.3624485
confounded$confounded_pair
#> [1] 3 4
```

``` r

parce_result <- lingam_parce(confounded$data, reg_method = "ols")
print(parce_result)
#> Bottom-Up ParceLiNGAM Result
#>   Variables : 6
#>   Independence measure: hsic
#>   Causal order: (x2, x3) -> x0 -> x4 -> x5 -> x1
#>   (NA entries in the adjacency matrix = unresolved order / suspected latent confounding)
#> 
#> Adjacency matrix (row = to, col = from):
#>       x0 x1     x2     x3    x4     x5
#> x0 0.000  0 -0.010  0.516 0.000  0.000
#> x1 0.479  0  0.447  0.060 0.025 -0.049
#> x2 0.000  0  0.000     NA 0.000  0.000
#> x3 0.000  0     NA  0.000 0.000  0.000
#> x4 0.497  0 -0.490 -0.001 0.000  0.000
#> x5 0.436  0  0.068  0.023 0.050  0.000
```

The causal order’s first element is the unresolved block, shown in
parentheses; here it correctly contains `x2` and `x3`. The corresponding
entries of the adjacency matrix are `NA`, while edges among the
remaining, fully-resolved variables are estimated as usual:

``` r

parce_result$causal_order[[1]]
#> [1] 3 4
parce_result$adjacency_matrix[confounded$confounded_pair, confounded$confounded_pair]
#>    x2 x3
#> x2  0 NA
#> x3 NA  0
```

Because a confounded variable’s true parents cannot be identified,
[`estimate_total_effect_parce()`](https://morimotoosamu.github.io/lingamr/reference/estimate_total_effect_parce.md)
warns and returns `NA` when asked for a total effect *from* a variable
in the unresolved block, but still computes normal estimates for
well-identified pairs:

``` r

# from a confounded variable: warns and returns NA
estimate_total_effect_parce(confounded$data, parce_result,
  from_index = confounded$confounded_pair[1], to_index = "x1"
)
#> Warning in estimate_total_effect_parce(confounded$data, parce_result,
#> from_index = confounded$confounded_pair[1], : x2 is part of an unresolved
#> causal order (suspected latent confounding); total effect cannot be estimated.
#> [1] NA

# a well-identified pair: a normal numeric estimate
estimate_total_effect_parce(confounded$data, parce_result,
  from_index = "x0", to_index = "x5"
)
#> [1] 0.5121874
```

[`get_error_independence_p_values_parce()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values_parce.md)
complements this by testing, for every pair of variables, whether their
residuals are independent — using the HSIC gamma-approximation test
([`hsic_test_gamma()`](https://morimotoosamu.github.io/lingamr/reference/hsic_test_gamma.md))
rather than a correlation test. Any pair touching a variable in the
unresolved block cannot be tested and comes back as `NA`;
well-identified pairs get a normal p-value:

``` r

round(get_error_independence_p_values_parce(confounded$data, parce_result), 3)
#>       x0    x1 x2 x3    x4    x5
#> x0    NA 0.453 NA NA 0.795 0.367
#> x1 0.453    NA NA NA 0.924 0.471
#> x2    NA    NA NA NA    NA    NA
#> x3    NA    NA NA NA    NA    NA
#> x4 0.795 0.924 NA NA    NA 0.873
#> x5 0.367 0.471 NA NA 0.873    NA
```

[`lingam_parce_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce_bootstrap.md)
provides bootstrap stability estimates in the same style as
[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md).
`NA` (unresolved) edges are treated as absent when aggregating, so
[`get_probabilities()`](https://morimotoosamu.github.io/lingamr/reference/get_probabilities.md)
and the other `BootstrapResult` query functions work as usual;
[`get_causal_order_stability()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_order_stability.md)
is the one exception, since ParceLiNGAM’s blocked causal order does not
fit its fixed-length format.

``` r

parce_bs <- lingam_parce_bootstrap(confounded$data,
  n_sampling = 10L, reg_method = "ols", seed = 1, verbose = FALSE
)
get_probabilities(parce_bs)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]  0.0  0.2  0.5  0.4  0.0  0.0
#> [2,]  0.5  0.0  0.5  0.4  0.2  0.3
#> [3,]  0.0  0.0  0.0  0.0  0.0  0.0
#> [4,]  0.1  0.2  0.2  0.0  0.1  0.0
#> [5,]  0.7  0.5  0.7  0.6  0.0  0.3
#> [6,]  0.6  0.4  0.6  0.6  0.4  0.0
```

## RCD (Repetitive Causal Discovery)

[`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)
(Repetitive Causal Discovery; Maeda and Shimizu 2020) tackles the same
latent-confounder problem as
[`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md),
but from a different angle: rather than searching for a causal order and
giving up on an **unresolved block** once a test is rejected, RCD
directly estimates each variable’s **ancestor set** and then checks
individual, parent-free pairs for a shared latent confounder. This makes
RCD’s output pair-level (which specific pairs are confounded) rather
than block-level (which set of variables could not be ordered).

[`generate_rcd_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_rcd_sample.md)
generates a 7-variable model in which `x6` is an unobserved common cause
of `x2` and `x4`; only `x0`-`x5` are returned as data.

``` r

# HSIC is O(n^2), so a moderate n keeps this article fast to build
rcd_confounded <- generate_rcd_sample(n = 300, seed = 1)
head(rcd_confounded$data)
#>            x0           x1         x2           x3         x4            x5
#> 1 -0.96839700 -0.023398020 -0.7514703 -0.473147026 -0.9674898 -0.0307310344
#> 2  1.31480945  0.424388510  1.0389688  0.001304296  1.4111680  0.0007741684
#> 3 -0.85973904 -0.025330472 -1.1731736 -0.034157640 -1.3389435 -0.0729373397
#> 4 -0.76463976  0.324413350 -0.7116352  0.078719765 -0.3333176  0.5074829194
#> 5  0.08152383  0.002364107 -0.2546673  0.036715357 -0.2040856  0.0044720536
#> 6 -0.30389184 -0.225029302 -0.5032249 -0.172267536  0.5221642 -0.0690391705
rcd_confounded$confounded_pair
#> [1] 3 5
```

``` r

rcd_result <- lingam_rcd(rcd_confounded$data)
print(rcd_result)
#> RCD Result
#>   Variables : 6
#> 
#> Ancestor sets:
#>   M(x0) = {x1, x3, x5}
#>   M(x1) = {x5}
#>   M(x2) = {x0, x1, x3, x5}
#>   M(x3) = {x5}
#>   M(x4) = {x0, x1, x3, x5}
#>   M(x5) = {}
#> 
#>   (NA entries in the adjacency matrix = suspected shared latent confounder)
#> 
#> Adjacency matrix (row = to, col = from):
#>       x0    x1 x2    x3 x4    x5
#> x0 0.000 1.116  0 0.989  0 0.000
#> x1 0.000 0.000  0 0.000  0 0.588
#> x2 0.810 0.000  0 0.000 NA 0.000
#> x3 0.000 0.000  0 0.000  0 0.449
#> x4 1.015 0.000 NA 0.000  0 0.000
#> x5 0.000 0.000  0 0.000  0 0.000
```

`ancestors_list` gives each variable’s estimated ancestors (not a causal
order), and the confounded pair’s adjacency-matrix entries are `NA`:

``` r

rcd_result$ancestors_list
#> $x0
#> [1] 2 4 6
#> 
#> $x1
#> [1] 6
#> 
#> $x2
#> [1] 1 2 4 6
#> 
#> $x3
#> [1] 6
#> 
#> $x4
#> [1] 1 2 4 6
#> 
#> $x5
#> integer(0)
rcd_result$adjacency_matrix[rcd_confounded$confounded_pair, rcd_confounded$confounded_pair]
#>    x2 x4
#> x2  0 NA
#> x4 NA  0
```

As with ParceLiNGAM,
[`estimate_total_effect_rcd()`](https://morimotoosamu.github.io/lingamr/reference/estimate_total_effect_rcd.md)
warns and returns `NA` when asked for a total effect *from* a confounded
variable:

``` r

# from a confounded variable: warns and returns NA
estimate_total_effect_rcd(rcd_confounded$data, rcd_result,
  from_index = rcd_confounded$confounded_pair[1], to_index = rcd_confounded$confounded_pair[2]
)
#> Warning in estimate_total_effect_rcd(rcd_confounded$data, rcd_result,
#> from_index = rcd_confounded$confounded_pair[1], : x2 is part of a suspected
#> latent confounder pair; total effect cannot be estimated.
#> [1] NA

# a well-identified pair: a normal numeric estimate
estimate_total_effect_rcd(rcd_confounded$data, rcd_result,
  from_index = "x5", to_index = "x0"
)
#>      x5 
#> 1.05674
```

[`get_error_independence_p_values_rcd()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values_rcd.md)
provides the same residual-independence check for RCD results: pairs
touching a confounded variable are `NA`, while the rest get a normal
p-value:

``` r

round(get_error_independence_p_values_rcd(rcd_confounded$data, rcd_result), 3)
#>       x0    x1 x2    x3 x4    x5
#> x0    NA 0.000 NA 0.494 NA 0.088
#> x1 0.000    NA NA 0.808 NA 0.401
#> x2    NA    NA NA    NA NA    NA
#> x3 0.494 0.808 NA    NA NA 0.017
#> x4    NA    NA NA    NA NA    NA
#> x5 0.088 0.401 NA 0.017 NA    NA
```

[`lingam_rcd_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd_bootstrap.md)
provides bootstrap stability estimates; the usual `BootstrapResult`
query functions apply.

## ParceLiNGAM or RCD?

| Aspect | ParceLiNGAM | RCD |
|----|----|----|
| Output granularity | Block-level: a set of variables that could not be ordered | Pair-level: which specific pairs are confounded |
| Main result | Blocked causal order + adjacency matrix with `NA` block | Ancestor sets + adjacency matrix with `NA` pairs |
| Question it answers | “Which part of the ordering can I trust?” | “Exactly which pairs share a hidden common cause?” |

Both rely on HSIC independence tests, whose cost grows as $`O(n^2)`$ —
with large $`n`$, consider subsampling. If relationships are nonlinear,
use
[CAM-UV](https://morimotoosamu.github.io/lingamr/articles/nonlinear.md)
instead.

## Related Articles

- [Method selection
  guide](https://morimotoosamu.github.io/lingamr/articles/method-selection.md)
  — which method fits your data
- [Nonlinear methods (RESIT and
  CAM-UV)](https://morimotoosamu.github.io/lingamr/articles/nonlinear.md)
  — CAM-UV handles latent confounders in nonlinear models
- [Bootstrap and
  diagnostics](https://morimotoosamu.github.io/lingamr/articles/bootstrap-diagnostics.md)
