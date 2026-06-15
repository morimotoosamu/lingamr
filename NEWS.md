# lingamr 0.1.0

* Initial CRAN submission.
* Direct LiNGAM (`lingam_direct()`) with selectable regression backends for
  adjacency-matrix estimation via `reg_method`: ordinary least squares
  (`"ols"`), LASSO (`"lasso"`), adaptive LASSO (`"adaptive_lasso"`), and ridge
  regression (`"ridge"`).
* `lingam_direct_bootstrap()` provides bootstrap stability assessment,
  including causal-order stability, and supports multi-core execution through
  the `parallel` and `n_cores` arguments (via `parallel::makePSOCKcluster()`).
  Sequential execution remains the default. Parallel runs use L'Ecuyer
  parallel RNG streams, so results are reproducible for a given
  `seed`/`n_cores` but differ numerically from the sequential path.
* Model diagnostics: residual independence and normality tests, plus a
  one-call `summary_lingam()`.
* Visualization with DiagrammeR (interactive) and ggplot2 `autoplot()`
  (static).
* broom-style tidiers (`tidy()` / `glance()`).
