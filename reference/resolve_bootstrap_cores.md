# Resolve the number of parallel cores for a bootstrap run

Shared by all `*_bootstrap()` functions. Applies the common policy:
`n_cores = NULL` defaults to a maximum of 2 cores for safety, the
resolved value is capped by the available cores and by `n_sampling`, and
a resolved value of 1 demotes the run to sequential execution.

## Usage

``` r
resolve_bootstrap_cores(parallel, n_cores, n_sampling)
```

## Arguments

- parallel:

  Whether parallel execution was requested (logical).

- n_cores:

  Requested number of cores (integer or NULL).

- n_sampling:

  Number of bootstrap iterations (validated positive integer).

## Value

`list(parallel = <logical>, n_cores = <integer or NULL>)` with the
possibly-demoted `parallel` flag and the resolved core count.
