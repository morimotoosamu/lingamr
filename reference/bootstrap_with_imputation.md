# Bootstrap with Multiple Imputation for Direct LiNGAM

Causal discovery on data containing missing values (NA). Each bootstrap
resample (drawn with replacement, missing values retained) is multiply
imputed into `n_repeats` complete datasets, and a common causal
structure is jointly estimated across those datasets with
[`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md)
(Shimizu 2012), treating the imputed copies as "groups" sharing one
causal order. R port of the Python
`lingam.tools.bootstrap_with_imputation()`.

## Usage

``` r
bootstrap_with_imputation(
  X,
  n_sampling,
  n_repeats = 10L,
  imputer = NULL,
  cd_fit = NULL,
  prior_knowledge = NULL,
  apply_prior_knowledge_softly = FALSE,
  seed = NULL,
  verbose = TRUE,
  parallel = FALSE,
  n_cores = NULL
)
```

## Arguments

- X:

  A numeric matrix or data frame (n_samples x n_features). May contain
  `NA`. If `X` has no missing values, a warning suggests using
  [`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md)
  instead, and estimation proceeds anyway.

- n_sampling:

  Number of bootstrap iterations (positive integer)

- n_repeats:

  Number of imputed datasets generated per bootstrap sample (positive
  integer, default `10L`). Ignored when a custom `imputer` is supplied;
  the number of datasets it returns is used instead.

- imputer:

  `NULL`, or a `function(X_boot)` returning a list of complete (no-NA)
  numeric matrices, each with the same dimensions as `X_boot`. Defaults
  to multiple imputation via `mice::mice(method = "norm")`.

- cd_fit:

  `NULL`, or a `function(X_list)` returning
  `list(causal_order = <integer vector, 1-based permutation>, adjacency_matrices = <list of p x p matrices, one per element of X_list>)`.
  Defaults to joint estimation via
  [`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md).

- prior_knowledge:

  Prior knowledge matrix (NULL allowed). Only used when `cd_fit = NULL`;
  a warning is issued if supplied together with a custom `cd_fit`.

- apply_prior_knowledge_softly:

  Apply prior knowledge softly (logical). Same restriction as
  `prior_knowledge`.

- seed:

  Random seed (NULL allowed). In sequential mode it is set once before
  the bootstrap loop and governs both the resampling and (via the global
  RNG) `mice`'s imputation; in parallel mode it seeds the L'Ecuyer
  parallel random-number streams (see Details).

- verbose:

  Whether to display progress (logical)

- parallel:

  Whether to use parallel processing (logical)

- n_cores:

  Number of cores to use (integer, NULL allowed). When `NULL`, the
  number of cores is limited to a maximum of 2 for safety. Ignored when
  `parallel = FALSE`.

## Value

An `ImputationBootstrapResult` (list) containing:

- `causal_orders`: `n_sampling` x `p` integer matrix (1-based causal
  order per iteration).

- `adjacency_matrices`: `array(n_sampling, n_repeats, p, p)`;
  `[, , i, j]` follows the lingamr convention (`B[i, j]` = coefficient
  of j -\> i).

- `resampled_indices`: `n_sampling` x `n` integer matrix of resampled
  row indices.

- `imputation_results`: `array(n_sampling, n_repeats, n, p)`; non-`NA`
  only at positions that were missing in that iteration's bootstrap
  resample.

## Details

**Procedure:** for each of `n_sampling` iterations, (1) resample `X`
with replacement (missing values are retained), (2) impute the resample
into `n_repeats` complete datasets, (3) jointly estimate one causal
structure shared by all `n_repeats` datasets with
[`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md).
This assumes the same causal structure underlies every imputed copy.

**Default imputer.** `mice::mice(method = "norm")` (Bayesian linear
regression, multiple imputation by chained equations) is the closest
standard R analogue of the upstream Python default
(`IterativeImputer(sample_posterior = TRUE)`). The two do not produce
numerically identical imputations.

**Custom `imputer` / `cd_fit`.** Supply your own imputation or
causal-discovery routine by passing a function with the signature
described above; the return value is validated and a descriptive error
is raised on violation. This replaces the abstract base classes
(`BaseMultipleImputation`, `BaseMultiGroupCDModel`) of the Python
original.

**Downstream analysis.** The result's shape (an extra `n_repeats`
dimension for `adjacency_matrices` and `imputation_results`) differs
from
[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md)'s
`BootstrapResult`, so it cannot be passed directly to
[`get_probabilities()`](https://morimotoosamu.github.io/lingamr/reference/get_probabilities.md)
etc. Use
[`as_bootstrap_result()`](https://morimotoosamu.github.io/lingamr/reference/as_bootstrap_result.md)
to collapse the `n_repeats` dimension (aggregating by median or mean)
into a `BootstrapResult`.

**On iteration failures:** each iteration is wrapped in
[`tryCatch()`](https://rdrr.io/r/base/conditions.html); a failing
iteration (e.g. `mice` fails to converge on a particular resample) is
skipped with a warning, and only if every iteration fails is an error
raised, mirroring
[`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md).

**On reproducibility:** During parallel execution, L'Ecuyer parallel
random number streams via
[`parallel::clusterSetRNGStream()`](https://rdrr.io/r/parallel/RngStream.html)
are used. Results are reproducible given the same `seed` and same
`n_cores`, but they do not numerically match the results of sequential
execution (`parallel = FALSE`). If you need results that exactly match
the sequential version, use `parallel = FALSE`. The 'mice' imputation
draws from each worker's global RNG, so it is covered by the same stream
mechanism.

## Examples

``` r
set.seed(1)
sample6 <- generate_lingam_sample_6(n = 300, seed = 1)
X <- sample6$data
X$x5[sample.int(nrow(X), size = round(0.1 * nrow(X)))] <- NA # MCAR 10% on x5

# \donttest{
if (requireNamespace("mice", quietly = TRUE)) {
  res <- bootstrap_with_imputation(X,
    n_sampling = 5L, n_repeats = 3L, seed = 42, verbose = FALSE
  )
  print(res)

  # Collapse the n_repeats dimension to reuse the existing bootstrap tooling
  bs <- as_bootstrap_result(res, aggregate = "median")
  get_probabilities(bs)
}
#> ImputationBootstrapResult: 5 samplings x 3 repeats, 6 features, 30 missing cells (original data)
#>     x0 x1 x2  x3  x4 x5
#> x0 0.0  0  0 0.8 0.2  0
#> x1 1.0  0  1 0.0 0.0  0
#> x2 0.0  0  0 1.0 0.0  0
#> x3 0.0  0  0 0.0 0.0  0
#> x4 0.8  0  1 0.2 0.0  0
#> x5 1.0  0  0 0.0 0.0  0
# }
```
