# Changelog

## lingamr (development version)

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
