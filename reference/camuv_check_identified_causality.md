# Check that no pair within a variable set already has a parent relation

Faithful port of `CAMUV._check_identified_causality()`.

## Usage

``` r
camuv_check_identified_causality(vars, P)
```

## Arguments

- vars:

  integer vector of variable indices

- P:

  current parent list

## Value

TRUE if no pair in `vars` has an identified causal relation yet
