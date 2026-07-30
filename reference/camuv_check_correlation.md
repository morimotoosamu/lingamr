# Check that every candidate parent is in the child's neighborhood

Faithful port of `CAMUV._check_correlation()`.

## Usage

``` r
camuv_check_correlation(child, parents, N)
```

## Arguments

- child:

  candidate child index

- parents:

  candidate parent indices

- N:

  neighborhood list from
  [`camuv_get_neighborhoods()`](https://morimotoosamu.github.io/lingamr/reference/camuv_get_neighborhoods.md)

## Value

TRUE if all parents are dependent on the child
