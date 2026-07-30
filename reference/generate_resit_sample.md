# Generate sample data from a nonlinear additive noise model (for RESIT)

Generates a 4-variable nonlinear SEM whose structure can only be
recovered by a nonlinear method such as
[`lingam_resit()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit.md);
linear LiNGAM variants are not expected to work on this data.

## Usage

``` r
generate_resit_sample(n = 300L, seed = NULL)
```

## Arguments

- n:

  number of samples (default: 300)

- seed:

  random seed (default: NULL, i.e. do not reset the RNG state)

## Value

list with three elements:

- `data`: data.frame of the 4 observed variables (`x0`-`x3`).

- `adjacency_matrix`: the true 4x4 adjacency matrix following the
  `m[to, from]` convention. Entries are 0/1 edge indicators (not
  coefficients), directly comparable to the output of
  [`lingam_resit()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit.md).

- `causal_order`: the true causal order (1-based column positions,
  source first).

## Details

The data-generating process (all error terms `e()` are
`runif(n, -0.5, 0.5)`):


    x0 ~ runif(-1, 1)
    x1 = 3.0 * x0^2 + e()
    x2 = 2.0 * tanh(x1) + 0.8 * x0^3 + e()
    x3 = 1.5 * sin(x2) + e()

## Examples

``` r
nonlinear <- generate_resit_sample(n = 300, seed = 1)
head(nonlinear$data)
#>           x0          x1           x2           x3
#> 1 -0.4689827  0.83354648  1.596482936  1.831404164
#> 2 -0.2557522 -0.20891458  0.003539709  0.272152307
#> 3  0.1457067  0.05628747 -0.237588016 -0.580258311
#> 4  0.8164156  1.96115504  2.607512659  0.451737357
#> 5 -0.5966361  0.94314057  1.779071173  1.193345462
#> 6  0.7967794  2.39567131  2.846533624 -0.001835138
nonlinear$adjacency_matrix
#>    x0 x1 x2 x3
#> x0  0  0  0  0
#> x1  1  0  0  0
#> x2  1  1  0  0
#> x3  0  0  1  0
```
