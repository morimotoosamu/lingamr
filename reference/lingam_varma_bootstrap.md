# Bootstrap for VARMA-LiNGAM

Evaluates the statistical reliability of the estimated time-series DAG
by resampling. Like
[`lingam_var_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_var_bootstrap.md),
this uses a **residual bootstrap**: the VARMA model is fitted once on
the original data, the residuals are resampled with replacement, and a
new series is rebuilt by the VARMA recursion before re-estimating
VARMA-LiNGAM on it. Port of the Python reference
`VARMALiNGAM.bootstrap`.

## Usage

``` r
lingam_varma_bootstrap(
  X,
  n_sampling,
  order = c(1L, 1L),
  criterion = "bic",
  measure = "pwling",
  reg_method = "adaptive_lasso",
  lambda = "BIC",
  init_method = "ols",
  prune = TRUE,
  seed = NULL,
  verbose = TRUE,
  parallel = FALSE,
  n_cores = NULL
)
```

## Arguments

- X:

  numeric matrix or data frame (n_samples x n_features), rows ordered in
  time.

- n_sampling:

  number of bootstrap iterations (positive integer).

- order:

  VARMA order `c(p, q)`. When `criterion` is not NULL, the order is
  selected once on the original data and then fixed across all
  iterations.

- criterion:

  order-selection criterion ("bic", "aic", "hqic") or NULL to use
  `order` directly.

- measure:

  independence measure for
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
  ("pwling"/"kernel").

- reg_method:

  regression method for the instantaneous matrix.

- lambda:

  penalty selection (see
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)).

- init_method:

  initial-weight method for adaptive LASSO.

- prune:

  logical; passed to
  [`lingam_varma()`](https://morimotoosamu.github.io/lingamr/reference/lingam_varma.md)
  on each iteration (default TRUE).

- seed:

  random seed (NULL allowed).

- verbose:

  whether to print progress (logical).

- parallel:

  whether to distribute iterations across cores (logical).

- n_cores:

  number of cores (integer or NULL; NULL caps at 2 for safety).

## Value

a `VARMABootstrapResult` object.

## Details

Reproducibility follows the same rules as
[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md):
with `parallel = TRUE`, L'Ecuyer streams via
[`parallel::clusterSetRNGStream()`](https://rdrr.io/r/parallel/RngStream.html)
make results reproducible for a given `seed` and `n_cores`, but they do
not match the sequential (`parallel = FALSE`) results.

**On iteration failures:** as in
[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md),
each iteration runs inside a
[`tryCatch()`](https://rdrr.io/r/base/conditions.html); a failing
iteration is reported as a warning and excluded from the result instead
of aborting the run. An error is raised only if every iteration fails.

As in the Python reference, the series regeneration omits the estimated
intercept, so resampled series are centered near zero even when the
original data are not; each refit re-estimates its own intercept, so the
resampled coefficient estimates are unaffected.

Total effects are estimated by the back-door regression of
[`estimate_varma_total_effect()`](https://morimotoosamu.github.io/lingamr/reference/estimate_varma_total_effect.md)
(the Python reference does the same) and cover the instantaneous block
and the AR lags 1..p; the MA (omega) blocks describe effects of past
disturbances, not of observed variables, and are therefore excluded from
`total_effects`.

## Examples

``` r
s <- generate_varmalingam_sample(n = 300, seed = 42)

# Fast example: OLS instantaneous structure, no pruning (no glmnet needed)
bs <- lingam_varma_bootstrap(s$data,
  n_sampling = 5L, order = c(1, 1), criterion = NULL,
  reg_method = "ols", prune = FALSE, seed = 1, verbose = FALSE
)
get_varma_probabilities(bs)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#> [1,]    0    0    0    1    1    1    1    1    1
#> [2,]    1    0    0    1    1    1    1    1    1
#> [3,]    1    1    0    1    1    1    1    1    1
```
