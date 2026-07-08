# Whether `xi` can be excluded from sink candidacy given the ancestor sets known so far

Faithful port of the exclusion rule in `rcd.py` (original 165-174
lines): `xi` cannot be the sink of `U` if (a) `xi` is already known to
be an ancestor of some other member of `U`, or (b) every other member of
`U` is already known to be an ancestor of `xi`.

## Usage

``` r
exists_ancestor_in_U(M, U, xi, xj_list)
```

## Arguments

- M:

  current ancestor-set list

- U:

  variable-set under consideration (unused except via `xj_list`, kept
  for signature parity with the upstream port)

- xi:

  candidate variable

- xj_list:

  `U` minus `xi`

## Value

TRUE if `xi` should be excluded from sink candidacy
