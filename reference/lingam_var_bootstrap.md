# Bootstrap for VAR-LiNGAM

Evaluates the statistical reliability of the estimated time-series DAG
by resampling. Unlike the i.i.d. row resampling used for Direct LiNGAM,
this uses a **residual bootstrap**: the VAR is fitted once on the
original data, the residuals are resampled with replacement, and a new
series is rebuilt by the VAR recursion before re-estimating VAR-LiNGAM
on it (this preserves the autoregressive structure). Port of the Python
reference `VARLiNGAM.bootstrap`.

## Usage

``` r
lingam_var_bootstrap(
  X,
  n_sampling,
  lags = 1L,
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

- lags:

  maximum lag order. When `criterion` is not NULL, the lag is selected
  once on the original data and then fixed across all iterations.

- criterion:

  lag-selection criterion ("bic", "aic", "hqic", "fpe") or NULL to use
  `lags` directly.

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
  [`lingam_var()`](https://morimotoosamu.github.io/lingamr/reference/lingam_var.md)
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

a `VARBootstrapResult` object.

## Details

Reproducibility follows the same rules as
[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md):
with `parallel = TRUE`, L'Ecuyer streams via
[`parallel::clusterSetRNGStream()`](https://rdrr.io/r/parallel/RngStream.html)
make results reproducible for a given `seed` and `n_cores`, but they do
not match the sequential (`parallel = FALSE`) results.

## Examples

``` r
s <- generate_varlingam_sample(n = 500, seed = 42)

# Fast example: OLS instantaneous structure, no pruning (no glmnet needed)
bs <- lingam_var_bootstrap(s$data,
  n_sampling = 10L, lags = 1, criterion = NULL,
  reg_method = "ols", prune = FALSE, seed = 1, verbose = FALSE
)
get_var_probabilities(bs)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    0    0    0    1    1    1
#> [2,]    1    0    0    1    1    1
#> [3,]    1    1    0    1    1    1
```
