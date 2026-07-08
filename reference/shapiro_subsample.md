# Deterministic subsample for the Shapiro-Wilk test

[`stats::shapiro.test()`](https://rdrr.io/r/stats/shapiro.test.html) has
a hard cap of 5000 observations, so larger inputs are thinned to
`SHAPIRO_MAX_N` evenly spaced values. The thinning is deterministic on
purpose: a random subsample would make the reported p-values change
between calls and silently consume the caller's RNG stream (breaking
downstream reproducibility). Even spacing over the input order does not
distort the marginal distribution being tested.

## Usage

``` r
shapiro_subsample(x)
```

## Arguments

- x:

  numeric vector

## Value

`x` itself when `length(x) <= SHAPIRO_MAX_N`, otherwise a deterministic
subsample of length `SHAPIRO_MAX_N`
