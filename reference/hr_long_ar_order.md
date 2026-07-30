# Choose the long autoregression order for Hannan-Rissanen

The first Hannan-Rissanen stage fits a long VAR(h) whose residuals proxy
the unobserved innovations. `h` follows the Box-Jenkins-style growth
rule `max(p + q, ceiling(1.5 (p + q)), floor(10 log10(n)))`, which
satisfies the consistency requirement that h grows with n while h^3/n
-\> 0.

## Usage

``` r
hr_long_ar_order(n, p_order, q_order, n_features)
```

## Arguments

- n:

  number of observations

- p_order:

  AR order

- q_order:

  MA order

- n_features:

  number of variables

## Value

the long AR order h (integer)
