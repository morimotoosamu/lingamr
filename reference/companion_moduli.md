# Eigenvalue moduli of the companion matrix

Builds the companion matrix of a lag-coefficient array (k, m, m) and
returns all eigenvalue moduli. Used for AR stationarity and MA
invertibility checks (both require every modulus to be below 1).

## Usage

``` r
companion_moduli(coefs)
```

## Arguments

- coefs:

  coefficient array (k, m, m)

## Value

numeric vector of eigenvalue moduli (empty when k = 0)
