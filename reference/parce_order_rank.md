# Find the position (rank) of a variable within a ParceLiNGAM causal order

All members of the unresolved block (if any) share the same rank (1).

## Usage

``` r
parce_order_rank(causal_order, idx)
```

## Arguments

- causal_order:

  list as produced by
  [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)

- idx:

  1-based variable index

## Value

integer rank, or NA if not found
