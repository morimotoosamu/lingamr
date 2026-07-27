# Report and drop failed bootstrap iterations

Shared post-processing for all `*_bootstrap()` functions: issues one
warning per failed iteration (an element with `ok = FALSE`), removes the
failures, and aborts only when every iteration failed.

## Usage

``` r
filter_bootstrap_failures(res_list)
```

## Arguments

- res_list:

  List of per-iteration results; each element has `ok` (logical) plus
  either the iteration payload or `iteration` / `message` describing the
  failure.

## Value

The list restricted to successful iterations (never empty).
