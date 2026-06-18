# Generate sample data from a VAR-LiNGAM model

Generates a 3-variable time series following a VAR-LiNGAM model with a
strictly acyclic instantaneous structure B0, a lag-1 coefficient matrix
M1, and non-Gaussian (uniform) errors.

## Usage

``` r
generate_varlingam_sample(n = 1000, seed = NULL)
```

## Arguments

- n:

  number of time points to return (after burn-in)

- seed:

  random seed (NULL allowed)

## Value

list with `data` (data frame, n x 3), `true_B0` (instantaneous matrix),
and `true_M1` (lag-1 coefficient matrix)

## Examples

``` r
sample <- generate_varlingam_sample(n = 500, seed = 1)
head(sample$data)
#>            x0          x1          x2
#> 1 -0.08134851  0.89990884 -1.13842754
#> 2 -0.66767312  0.67124732 -1.58403973
#> 3 -1.20176088 -0.01791795 -1.37283242
#> 4  0.09281405  0.40345059  0.06231193
#> 5  0.32280573  0.61074550  0.65849296
#> 6 -0.24691355 -0.81843074  1.25860250
```
