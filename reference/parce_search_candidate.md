# Candidate sink variables under prior knowledge (ParceLiNGAM direction)

Unlike
[`search_candidate()`](https://morimotoosamu.github.io/lingamr/reference/search_candidate.md)
(used by top-down DirectLiNGAM), ParceLiNGAM searches from the sink
side, so the filter is simply "exclude variables that appear as the
'from' side of a known partial order" (they cannot be a sink because
something is known to cause a variable through them... more precisely: a
variable known to *cause* another remaining variable cannot itself be
the next sink).

## Usage

``` r
parce_search_candidate(U, partial_orders)
```

## Arguments

- U:

  set of currently undetermined variables

- partial_orders:

  matrix of (from, to) pairs from
  [`extract_partial_orders()`](https://morimotoosamu.github.io/lingamr/reference/extract_partial_orders.md)

## Value

candidate index vector (subset of U, or U itself if no candidates
remain)
