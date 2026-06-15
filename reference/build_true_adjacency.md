# Build the true adjacency matrix from edge specifications

Build the true adjacency matrix from edge specifications

## Usage

``` r
build_true_adjacency(var_names, from, to, coef)
```

## Arguments

- var_names:

  vector of variable names

- from:

  vector of cause (source) variable names for the edges

- to:

  vector of effect (target) variable names for the edges (same length as
  `from`)

- coef:

  vector of edge coefficients (same length as `from`)

## Value

adjacency matrix (p x p). `m[to, from] = coef` (row = to, col = from)
