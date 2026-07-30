# Changelog

## lingamr (development version)

- Added
  [`lingam_camuv()`](https://morimotoosamu.github.io/lingamr/reference/lingam_camuv.md)
  and
  [`generate_camuv_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_camuv_sample.md),
  an R port of CAM-UV (Causal Additive Models with Unobserved Variables;
  Maeda and Shimizu
  2021. for causal discovery on nonlinear additive models with
        unobserved variables. The result is a new `CAMUVResult` class:
        the adjacency matrix holds 0/1 edge indicators (no coefficients
        are estimated), and variable pairs suspected to be connected
        through an unobserved causal or backdoor path are `NA`, with
        matching [`print()`](https://rdrr.io/r/base/print.html) /
        [`tidy()`](https://generics.r-lib.org/reference/tidy.html) /
        [`glance()`](https://generics.r-lib.org/reference/glance.html) /
        [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
        methods. The regressor is pluggable as in
        [`lingam_resit()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit.md)
        (the Python implementation hardcodes pygam’s `LinearGAM`; the
        default here is mgcv’s `gam()`). Prior knowledge uses the
        upstream pair format (`c(i, j)` = “variable i cannot be a cause
        of variable j”, 1-based). `independence = "fcorr"` is restricted
        to `num_explanatory_vals = 2`, where the Python implementation
        silently breaks on larger subsets. There is no bootstrap
        variant, matching the Python implementation.
- Added
  [`lingam_resit()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit.md),
  [`lingam_resit_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit_bootstrap.md),
  and
  [`generate_resit_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_resit_sample.md),
  an R port of RESIT (Regression with Subsequent Independence Test;
  Peters et al. 2014) for causal discovery on nonlinear additive noise
  models. The regressor is pluggable: the default `"gam"` fits
  smoothing-spline GAMs via the suggested package mgcv, and any function
  `function(X, y)` returning fitted values can be supplied instead. The
  result is a new `ResitResult` class whose adjacency matrix holds 0/1
  edge indicators (RESIT estimates no coefficients), with matching
  [`print()`](https://rdrr.io/r/base/print.html) /
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) /
  [`glance()`](https://generics.r-lib.org/reference/glance.html) /
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  methods. Unlike the Python implementation, the bootstrap stores no
  all-zero total-effects array (`total_effects` is `NULL`); total
  effects are undefined for nonlinear models.
- [`hsic_test_gamma()`](https://morimotoosamu.github.io/lingamr/reference/hsic_test_gamma.md)
  (internal) now accepts matrix arguments, combining the columns into a
  single multivariate Gaussian kernel as in the Python implementation;
  univariate callers are unchanged (bit-identical results).
- Fixed
  [`setup_cluster_worker()`](https://morimotoosamu.github.io/lingamr/reference/setup_cluster_worker.md)
  (internal, shared by all `*_bootstrap()` functions) to fall back to
  exporting the development namespace when the installed copy of the
  package is stale, instead of silently mixing old and new code on the
  workers.
- [`lingam_lim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_lim.md)
  gains an `is_poisson` argument for Poisson-type count discrete
  variables. With `is_poisson = TRUE`, discrete columns
  (`is_continuous = FALSE`) are validated as non-negative integer counts
  and the local search phase scores them with Poisson regression
  log-likelihoods (unregularized multivariate `glm(family = poisson())`
  with an intercept; closed-form intercept-only Poisson MLE for
  parentless variables). This is an intent-faithful port of the Python
  implementation’s `fit(is_poisson=True)`: the global optimization phase
  keeps the logistic surrogate loss, and several upstream scoring bugs
  are deliberately not reproduced (see the Details section of
  [`?lingam_lim`](https://morimotoosamu.github.io/lingamr/reference/lingam_lim.md)).
  [`generate_lim_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_lim_sample.md)
  gains a matching `is_poisson` argument, and `LiMResult` objects now
  carry an `is_poisson` field.

## lingamr 0.1.2

CRAN release: 2026-07-17

- Wrapped long-running bootstrap examples in `\donttest{}` to avoid CRAN
  NOTE for elapsed time \> 10s (`lingam_parce_bootstrap`,
  `lingam_rcd_bootstrap`).

## lingamr 0.1.1

- Extended the broom tidiers and
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  to the new result classes:
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) methods for
  `LiMResult`, `ParceLingamResult` and `RCDResult` (keeping `NA`
  adjacency entries visible as `estimate = NA` rows),
  `MultiGroupLingamResult` and `MultiGroupBootstrapResult` (stacked with
  a `group` column), and `ImputationBootstrapResult` (collapsed via
  [`as_bootstrap_result()`](https://morimotoosamu.github.io/lingamr/reference/as_bootstrap_result.md));
  [`glance()`](https://generics.r-lib.org/reference/glance.html) methods
  for `LiMResult` (`n_discrete`), `ParceLingamResult` (`n_na_entries`),
  `RCDResult` (`n_confounded_pairs`), and `MultiGroupLingamResult`
  (`n_groups`); and
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  methods for `LiMResult`, `ParceLingamResult`, `RCDResult` (suspected
  latent-confounder / unresolved pairs drawn as dashed segments), and
  `MultiGroupLingamResult` (one group at a time via the `group`
  argument).
- Added
  [`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md),
  [`lingam_rcd_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd_bootstrap.md),
  and
  [`generate_rcd_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_rcd_sample.md),
  an R port of RCD (Repetitive Causal Discovery; Maeda and Shimizu 2020)
  for causal discovery robust against latent confounders. Unlike
  [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md),
  RCD does not recover a causal order; instead it repeatedly extracts
  each variable’s ancestor set (`ancestors_list`), narrows ancestor sets
  down to direct parents, and tests remaining parent-free pairs for a
  shared latent confounder, marking the corresponding adjacency-matrix
  entries `NA`.
  [`estimate_total_effect_rcd()`](https://morimotoosamu.github.io/lingamr/reference/estimate_total_effect_rcd.md)
  and
  [`get_error_independence_p_values_rcd()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values_rcd.md)
  are the `RCDResult` counterparts of
  [`estimate_total_effect()`](https://morimotoosamu.github.io/lingamr/reference/estimate_total_effect.md)
  and
  [`get_error_independence_p_values()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values.md).
  Reuses the HSIC and F-correlation independence measures added for
  [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md),
  and adds an optional `MLHSICR` regression mode (HSIC-sum minimization
  via `stats::optim(method = "L-BFGS-B")`) as a fallback when OLS
  residuals are not independent of the explanatory variables.
- Added
  [`bootstrap_with_imputation()`](https://morimotoosamu.github.io/lingamr/reference/bootstrap_with_imputation.md),
  an R port of the Python `lingam.tools.bootstrap_with_imputation()`,
  for causal discovery on data containing missing values. Each bootstrap
  resample (drawn with replacement, missing values retained) is multiply
  imputed into `n_repeats` complete datasets (by default via
  `mice::mice(method = "norm")`, a new `Suggests` dependency), and a
  common causal structure shared by all imputed datasets is jointly
  estimated with
  [`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md).
  Imputation and causal-discovery estimation can be swapped for custom
  implementations via the `imputer` and `cd_fit` arguments; their return
  values are validated with descriptive errors on violation. The result
  is an `ImputationBootstrapResult`, whose extra `n_repeats` dimension
  can be collapsed into a regular `BootstrapResult` with the new
  [`as_bootstrap_result()`](https://morimotoosamu.github.io/lingamr/reference/as_bootstrap_result.md)
  helper to reuse the existing bootstrap analysis functions
  ([`get_probabilities()`](https://morimotoosamu.github.io/lingamr/reference/get_probabilities.md),
  [`get_causal_direction_counts()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_direction_counts.md),
  etc.).
- Added
  [`evaluate_model_fit()`](https://morimotoosamu.github.io/lingamr/reference/evaluate_model_fit.md),
  an R port of the Python `lingam.utils.evaluate_model_fit()`. Fits the
  causal graph implied by an estimated adjacency matrix (or a lingamr
  result object such as `LingamResult` / `ParceLingamResult` /
  `LiMResult`) as a structural equation model via
  [`lavaan::sem()`](https://rdrr.io/pkg/lavaan/man/sem.html) (a new
  `Suggests` dependency) and returns standard SEM fit measures (CFI,
  RMSEA, AIC/BIC, etc.). `NA` entries marking a latent confounder pair
  are represented as a residual covariance in the lavaan model,
  equivalent to the latent-variable representation used by the Python
  `semopy`-based original.
- Added
  [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md),
  [`lingam_parce_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce_bootstrap.md),
  and
  [`generate_parce_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_parce_sample.md),
  an R port of BottomUpParceLiNGAM (Tashiro et al. 2014) for causal
  discovery robust against latent confounders. The algorithm searches
  for a causal order from the sink side and stops once an independence
  test is rejected; variables it could not order are returned as a
  single unresolved block, and the corresponding adjacency-matrix
  entries are `NA`.
  [`estimate_total_effect_parce()`](https://morimotoosamu.github.io/lingamr/reference/estimate_total_effect_parce.md)
  and
  [`get_error_independence_p_values_parce()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values_parce.md)
  are the `ParceLingamResult` counterparts of
  [`estimate_total_effect()`](https://morimotoosamu.github.io/lingamr/reference/estimate_total_effect.md)
  and
  [`get_error_independence_p_values()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values.md).
  Adds two new internal-only independence measures reusable by future
  ports: an HSIC gamma-approximation test (`R/hsic.r`) and F-correlation
  / kernel canonical correlation (`R/f_correlation.r`).
- Added
  [`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md),
  [`lingam_multi_group_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group_bootstrap.md),
  [`get_group_result()`](https://morimotoosamu.github.io/lingamr/reference/get_group_result.md),
  and
  [`generate_multi_group_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_multi_group_sample.md),
  an R port of MultiGroupDirectLiNGAM (Shimizu 2012) for jointly
  estimating a Direct LiNGAM model across multiple datasets (“groups”)
  that share a common causal order but may have different structural
  coefficients. Per-group analysis (total causal effects, independence
  tests, plotting) reuses the existing single-group functions via
  [`get_group_result()`](https://morimotoosamu.github.io/lingamr/reference/get_group_result.md),
  which extracts a group as a plain `LingamResult`.
- Added
  [`lingam_high_dim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_high_dim.md),
  an R port of HighDimDirectLiNGAM (Wang & Drton 2020) for causal
  discovery on high-dimensional data (large `p`, or `p > n`). Causal
  order search uses moment statistics of non-Gaussianity instead of
  pairwise independence measures, and is deterministic.
- Added
  [`lingam_lim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_lim.md)
  and
  [`generate_lim_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_lim_sample.md),
  an R port of the LiM (LiNGAM for Mixed data) algorithm (Zeng et
  al. 2022) for causal discovery on data containing a mixture of
  continuous and binary (0/1) discrete variables.
- Fixed a condition in the kernel-based independence measure
  (`measure = "kernel"`) that made soft prior knowledge silently
  ineffective.
- Fixed `reg_method = "ridge"` erroring inside
  [`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md)
  and
  [`estimate_total_effect()`](https://morimotoosamu.github.io/lingamr/reference/estimate_total_effect.md)
  /
  [`estimate_all_total_effects()`](https://morimotoosamu.github.io/lingamr/reference/estimate_all_total_effects.md).
- Fixed `lambda = "oracle"` not being rejected upfront for
  `reg_method = "lasso"` (only `"ridge"` was previously validated),
  which previously surfaced as an unclear `glmnet` error.
- Fixed a data-scale dependence in the default adaptive-LASSO
  regularization path (`fit_regression.r`): the AIC/BIC lambda search
  grid is now scaled to the response’s magnitude instead of using a
  fixed absolute grid.
- [`select_var_lag()`](https://morimotoosamu.github.io/lingamr/reference/select_var_lag.md)
  now guards against selecting an overfit, near-saturated lag order when
  the sample size is small relative to the number of variables and
  candidate lags.
- [`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md)
  no longer aborts entirely when a single bootstrap iteration fails
  (e.g. a degenerate resample); the failing iteration is now skipped
  with a warning, and results reflect however many iterations succeeded.
- Added a `compute_total_effects` argument to
  [`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md)
  to skip the (comparatively expensive) total-effects estimation step
  when only edge/order stability is needed.
- [`get_causal_direction_counts()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_direction_counts.md)
  is now vectorized and substantially faster for large bootstrap
  results.
- `get_error_independence_p_values(method = "kendall")` now warns for
  large `n`, where Kendall’s tau is O(n^2) per variable pair.
- The kernel-based independence measure (`measure = "kernel"`) now
  switches to an incomplete-Cholesky low-rank approximation for
  `n > 1000`, cutting per-pair cost from O(n^3) to about O(n\*d^2)
  (~200x faster at n = 5000); `n <= 1000` still uses the exact
  computation.
- Removed unconditional Suggests-package dependencies from examples, and
  added `\examples` to the remaining exported `print.*` methods.
- Expanded test coverage (previously untested `BootstrapResult` query
  functions, numerical validation of total-effect estimates, and
  additional input-validation tests).

## lingamr 0.1.0

- Initial CRAN submission.
- Direct LiNGAM
  ([`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md))
  with selectable regression backends for adjacency-matrix estimation
  via `reg_method`: ordinary least squares (`"ols"`), LASSO (`"lasso"`),
  adaptive LASSO (`"adaptive_lasso"`), and ridge regression (`"ridge"`).
- [`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md)
  provides bootstrap stability assessment, including causal-order
  stability, and supports multi-core execution through the `parallel`
  and `n_cores` arguments (via
  [`parallel::makePSOCKcluster()`](https://rdrr.io/r/parallel/makeCluster.html)).
  Sequential execution remains the default. Parallel runs use L’Ecuyer
  parallel RNG streams, so results are reproducible for a given
  `seed`/`n_cores` but differ numerically from the sequential path.
- Model diagnostics: residual independence and normality tests, plus a
  one-call
  [`summary_lingam()`](https://morimotoosamu.github.io/lingamr/reference/summary_lingam.md).
- Visualization with DiagrammeR (interactive) and ggplot2
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  (static).
- broom-style tidiers
  ([`tidy()`](https://generics.r-lib.org/reference/tidy.html) /
  [`glance()`](https://generics.r-lib.org/reference/glance.html)).
