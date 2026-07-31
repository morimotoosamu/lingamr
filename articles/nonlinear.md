# Nonlinear Methods: RESIT and CAM-UV

LiNGAM assumes that causal relationships are **linear**. When a cause
acts on its effect through a nonlinear function (saturation, thresholds,
quadratic dose-response, …), linear methods can pick the wrong direction
while giving no warning. This article covers the two nonlinear methods
in `lingamr`:

- **RESIT**
  ([`lingam_resit()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit.md)):
  nonlinear additive noise models, assuming *no* latent confounder.
- **CAM-UV**
  ([`lingam_camuv()`](https://morimotoosamu.github.io/lingamr/reference/lingam_camuv.md)):
  causal additive models that additionally allow **unobserved
  variables**.

Both methods return 0/1 edge indicators rather than coefficients — with
nonlinear causal functions, a single number cannot summarize “the effect
of $`x_j`$ on $`x_i`$”, so total-effect estimation is intentionally
unavailable.

Both also require the suggested package `mgcv` (for GAM regressions).

``` r

library(lingamr)
```

## When Linear LiNGAM Fails on Nonlinear Data

[`generate_resit_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_resit_sample.md)
generates a 4-variable chain whose causal functions are all nonlinear
($`x_1 = 3x_0^2 + e`$, $`x_2 = 2\tanh(x_1) + 0.8x_0^3 + e`$,
$`x_3 = 1.5\sin(x_2) + e`$):

``` r

nonlin <- generate_resit_sample(n = 300, seed = 1)
head(nonlin$data)
#>           x0          x1           x2           x3
#> 1 -0.4689827  0.83354648  1.596482936  1.831404164
#> 2 -0.2557522 -0.20891458  0.003539709  0.272152307
#> 3  0.1457067  0.05628747 -0.237588016 -0.580258311
#> 4  0.8164156  1.96115504  2.607512659  0.451737357
#> 5 -0.5966361  0.94314057  1.779071173  1.193345462
#> 6  0.7967794  2.39567131  2.846533624 -0.001835138

# True 0/1 adjacency matrix (row = to, column = from)
nonlin$adjacency_matrix
#>    x0 x1 x2 x3
#> x0  0  0  0  0
#> x1  1  0  0  0
#> x2  1  1  0  0
#> x3  0  0  1  0
```

Applying **linear** Direct LiNGAM to this data gives a poor result — the
estimated order need not match the true order `x0, x1, x2, x3`:

``` r

linear_fit <- lingam_direct(nonlin$data)
colnames(nonlin$data)[linear_fit$causal_order]
#> [1] "x3" "x2" "x1" "x0"
```

This is not a defect of Direct LiNGAM: its linearity assumption simply
does not hold here. A nonlinear method is required.

## RESIT: Nonlinear Additive Noise Models

**RESIT** (Regression with Subsequent Independence Test; Peters et
al. 2014) assumes a nonlinear additive noise model

``` math
x_i = f_i(\mathrm{parents}(x_i)) + e_i
```

and recovers the structure in two phases:

1.  **Order search**: repeatedly detach the most sink-like variable —
    the one whose regression residual on all remaining variables is
    least dependent on them (smallest HSIC statistic).
2.  **Pruning**: for each candidate parent, test whether the residual
    regressed on the other parents is already independent of the parent
    set; if so (HSIC p-value above `alpha`), drop the parent.

``` r

resit_result <- lingam_resit(nonlin$data)
print(resit_result)
#> RESIT Result
#>   Variables : 4
#>   Regressor : gam
#>   Causal order: x0 -> x1 -> x2 -> x3
#> 
#> Adjacency matrix (row = to, col = from):
#>   (entries are 0/1 edge indicators, not coefficients)
#>    x0 x1 x2 x3
#> x0  0  0  0  0
#> x1  1  0  0  0
#> x2  1  1  0  0
#> x3  0  0  1  0
```

The estimated causal order and the 0/1 adjacency matrix recover the true
nonlinear structure:

``` r

colnames(nonlin$data)[resit_result$causal_order]
#> [1] "x0" "x1" "x2" "x3"
resit_result$adjacency_matrix
#>    x0 x1 x2 x3
#> x0  0  0  0  0
#> x1  1  0  0  0
#> x2  1  1  0  0
#> x3  0  0  1  0
```

Since the entries are edge indicators,
[`plot_adjacency()`](https://morimotoosamu.github.io/lingamr/reference/plot_adjacency.md)
draws an unweighted graph:

``` r

resit_result$adjacency_matrix |>
  plot_adjacency(
    labels = colnames(nonlin$data),
    title  = "RESIT estimate (edges are 0/1, not coefficients)"
  )
```

### Choosing `alpha`

Note the direction of `alpha` (default 0.01): a parent is removed when
the HSIC p-value *exceeds* `alpha`, so **larger** values make the test
stricter about declaring independence and therefore **keep more edges**.

### Custom Regressors

All internal regressions use a smoothing-spline GAM (`mgcv`) by default.
Any nonlinear regressor can be substituted by passing a function
`function(X, y)` that returns fitted values — for example, random
forests:

``` r

rf_regressor <- function(X, y) {
  fit <- randomForest::randomForest(X, y)
  as.numeric(predict(fit, X))
}
lingam_resit(nonlin$data, regressor = rf_regressor)
```

### Bootstrap

[`lingam_resit_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit_bootstrap.md)
works like the other bootstrap functions; the aggregated “probabilities”
are detection frequencies of each 0/1 edge. RESIT performs $`O(p^2)`$
GAM fits and HSIC tests ($`O(n^2)`$ each) *per resample*, so keep
`n_sampling` moderate:

``` r

resit_bs <- lingam_resit_bootstrap(nonlin$data,
  n_sampling = 20L, seed = 1, verbose = FALSE
)
get_probabilities(resit_bs)
#>      [,1] [,2] [,3] [,4]
#> [1,] 0.00 0.00  0.0 0.00
#> [2,] 1.00 0.00  0.1 0.05
#> [3,] 1.00 0.90  0.0 0.00
#> [4,] 0.05 0.05  1.0 0.00
```

## CAM-UV: Nonlinear Models with Unobserved Variables

RESIT assumes all relevant variables are observed. **CAM-UV** (Maeda and
Shimizu 2021) drops that assumption: it identifies each variable’s
direct parents by scanning variable subsets and testing residual
independence, and reports variable pairs whose dependence cannot be
explained by any identified edge as connected through an unobserved
variable:

- **UCP** (unobserved causal path): a directed path through an
  unobserved intermediate variable.
- **UBP** (unobserved backdoor path): a common unobserved ancestor (a
  latent confounder).

Such pairs are marked `NA` in the adjacency matrix rather than being
forced into a possibly wrong orientation.

[`generate_camuv_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_camuv_sample.md)
follows the Python tutorial’s structure: among six observed variables,
`u1` (latent) confounds `x3` and `x4` (UBP), and `u2` (latent) mediates
the path from `x2` to `x5` (UCP):

``` r

camuv_dat <- generate_camuv_sample(n = 500, seed = 1)
head(camuv_dat$data)
#>           x0          x1         x2         x3        x4        x5
#> 1 -0.8276528  0.13274564  0.1053562 2.06976237 2.2118223 0.7967904
#> 2 -0.4513472  0.48681681  0.6321648 1.68111156 1.7569962 0.4177104
#> 3  0.2571408  0.94440653 -0.3991332 0.04983659 0.6192421 1.3267213
#> 4  1.4407966  2.52524438  1.5559126 0.33134200 1.9471524 0.3952331
#> 5 -1.0529335 -0.03798973 -1.3050977 2.38726402 0.9939968 1.8703122
#> 6  1.4061430  2.96254650 -1.5761295 1.77541198 1.4811392 3.7850018

# The UBP (x3-x4) and UCP (x2-x5) entries of the true matrix are NA
camuv_dat$adjacency_matrix
#>    x0 x1 x2 x3 x4 x5
#> x0  0  0  0  0  0  0
#> x1  1  0  0  0  0  0
#> x2  0  0  0  0  0 NA
#> x3  1  0  0  0 NA  0
#> x4  0  0  1 NA  0  0
#> x5  0  0 NA  0  0  0
camuv_dat$confounded_pairs
#>      var1 var2
#> [1,]    3    6
#> [2,]    4    5
```

``` r

camuv_result <- lingam_camuv(camuv_dat$data)
print(camuv_result)
#> CAM-UV Result
#>   Variables : 6
#>   Regressor : gam
#> 
#> Parent sets:
#>   P(x0) = {}
#>   P(x1) = {x0}
#>   P(x2) = {}
#>   P(x3) = {x0}
#>   P(x4) = {x2}
#>   P(x5) = {}
#> 
#> Pairs with an unobserved causal/backdoor path (UCP/UBP):
#>   x2 -- x5
#>   x3 -- x4
#> 
#> Adjacency matrix (row = to, col = from):
#>   (entries are 0/1 edge indicators, not coefficients;
#>    NA = pair connected through an unobserved variable)
#>    x0 x1 x2 x3 x4 x5
#> x0  0  0  0  0  0  0
#> x1  1  0  0  0  0  0
#> x2  0  0  0  0  0 NA
#> x3  1  0  0  0 NA  0
#> x4  0  0  1 NA  0  0
#> x5  0  0 NA  0  0  0
```

Like
[`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md),
there is no causal order; the result holds each variable’s direct
parents and the flagged pairs:

``` r

camuv_result$parents_list
#> $x0
#> integer(0)
#> 
#> $x1
#> [1] 1
#> 
#> $x2
#> integer(0)
#> 
#> $x3
#> [1] 1
#> 
#> $x4
#> [1] 3
#> 
#> $x5
#> integer(0)
camuv_result$confounded_pairs
#>      var1 var2
#> [1,]    3    6
#> [2,]    4    5
camuv_result$adjacency_matrix
#>    x0 x1 x2 x3 x4 x5
#> x0  0  0  0  0  0  0
#> x1  1  0  0  0  0  0
#> x2  0  0  0  0  0 NA
#> x3  1  0  0  0 NA  0
#> x4  0  0  1 NA  0  0
#> x5  0  0 NA  0  0  0
```

### Tuning

- `num_explanatory_vals` (default 2) bounds the size of the variable
  subsets scanned for parents. Larger values increase statistical power
  but the cost grows combinatorially.
- `independence = "fcorr"` switches the independence measure from the
  HSIC test to the F-correlation with threshold `ind_corr`; it is only
  supported with `num_explanatory_vals = 2`.
- Prior knowledge is given as pairs `c(i, j)` meaning “variable i cannot
  cause variable j” (1-based indices).

There is no bootstrap variant, matching the Python implementation.

## RESIT or CAM-UV?

| Aspect | RESIT | CAM-UV |
|----|----|----|
| Latent confounders | Not allowed | Detected and reported (UCP/UBP pairs) |
| Output | Causal order + 0/1 adjacency matrix | Parents list + adjacency matrix with `NA` pairs |
| Bootstrap | [`lingam_resit_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit_bootstrap.md) | Not available |
| Cost | $`O(p^2)`$ GAM fits + HSIC tests | Combinatorial in `num_explanatory_vals` |

Start with RESIT when you can argue that all major common causes are
observed; use CAM-UV when latent confounding is plausible. Both build
$`n \times n`$ Gram matrices per HSIC test, so neither is recommended
for $`n`$ in the thousands — subsample first if needed.

If the relationships are linear and you suspect latent confounders, the
linear methods in the [latent confounders
article](https://morimotoosamu.github.io/lingamr/articles/latent-confounders.md)
(ParceLiNGAM, RCD) are cheaper and return coefficient estimates.

## Related Articles

- [Method selection
  guide](https://morimotoosamu.github.io/lingamr/articles/method-selection.md)
  — which method fits your data
- [Latent confounders (ParceLiNGAM and
  RCD)](https://morimotoosamu.github.io/lingamr/articles/latent-confounders.md)
  — linear counterparts
- [Bootstrap and
  diagnostics](https://morimotoosamu.github.io/lingamr/articles/bootstrap-diagnostics.md)
