# Validate a prior knowledge matrix and convert unknowns to NA

Checks the shape and that every entry is -1 (unknown), 0 (no path), 1
(path), or NA (treated as unknown). Anything else (e.g. 0.5 or 2) would
otherwise be silently interpreted as "path exists" by the candidate
search, so it is rejected here.

## Usage

``` r
validate_prior_knowledge(prior_knowledge, n_features)
```

## Arguments

- prior_knowledge:

  prior knowledge matrix

- n_features:

  expected number of variables

## Value

the validated matrix with negative entries replaced by NA
