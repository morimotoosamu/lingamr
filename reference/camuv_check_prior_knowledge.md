# Check whether prior knowledge forbids a candidate parent set

Faithful port of `CAMUV._check_prior_knowledge()`. Returns TRUE when the
combination is *blocked* (some candidate parent is forbidden from being
a cause of the child).

## Usage

``` r
camuv_check_prior_knowledge(pk_forbidden, parents, child)
```

## Arguments

- pk_forbidden:

  forbidden-cause list from
  [`camuv_make_pk_dict()`](https://morimotoosamu.github.io/lingamr/reference/camuv_make_pk_dict.md),
  or NULL

- parents:

  candidate parent indices

- child:

  candidate child index

## Value

TRUE if blocked by prior knowledge
