# Join psi and omega arrays into a single wide matrix

Returns `cbind(psi_0, ..., psi_p, omega_1, ..., omega_q)` with shape (m
x m(1 + p + q)), keeping the `[i, j]` = j -\> i convention in each
block.

## Usage

``` r
joined_varma_matrix(psis, omegas)
```

## Arguments

- psis:

  array (1 + p, m, m)

- omegas:

  array (q, m, m)

## Value

matrix (m x m(1 + p + q))
