# Build the CAM-UV adjacency matrix

Faithful port of `CAMUV._estimate_adjacency_matrix()`:
`B[child, parent] = 1` (edge indicators, not coefficients, since the
causal functions are nonlinear), and both entries of each confounded
pair are `NA`.

## Usage

``` r
camuv_build_adjacency_matrix(X, P, U)
```

## Arguments

- X:

  data matrix (for dimnames)

- P:

  parent list from
  [`camuv_find_parents()`](https://morimotoosamu.github.io/lingamr/reference/camuv_find_parents.md)

- U:

  list of confounded pairs (length-2 integer vectors)

## Value

adjacency matrix B (n_features x n_features)
