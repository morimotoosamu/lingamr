# lingamr (development version)

# lingamr 0.1.0.9005

* `lingam_direct_bootstrap()` gained `parallel` and `n_cores` arguments for
  multi-core execution via `parallel::makePSOCKcluster()`. Sequential
  execution remains the default. Parallel runs use L'Ecuyer parallel RNG
  streams, so results are reproducible for a given `seed`/`n_cores` but differ
  numerically from the sequential path.
* Initial CRAN submission.
