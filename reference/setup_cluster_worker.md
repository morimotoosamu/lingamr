# Make this package's functions available on parallel cluster workers

Used by all `*_bootstrap()` functions before dispatching iterations via
[`parallel::parLapply()`](https://rdrr.io/r/parallel/clusterApply.html).
Tries to have each worker
[`library()`](https://rdrr.io/r/base/library.html) the installed
package; if that is not possible (e.g. during `devtools::load_all()`
development, where the package is not installed), falls back to
exporting every object in the namespace environment of a representative
function from the algorithm being bootstrapped.

## Usage

``` r
setup_cluster_worker(cl, fun)
```

## Arguments

- cl:

  A `parallel` cluster object (from
  [`parallel::makePSOCKcluster()`](https://rdrr.io/r/parallel/makeCluster.html)).

- fun:

  A function belonging to the package whose namespace environment should
  be used as the fallback export source (e.g. `lingam_direct`,
  `lingam_rcd`, `lingam_parce`, `lingam_multi_group`).

## Value

`NULL`, invisibly. Called for the side effect of preparing `cl`.
