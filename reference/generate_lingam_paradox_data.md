# Generate Paradoxical Data Where DirectLiNGAM Struggles

Generates a synthetic dataset designed to favor ICA-LiNGAM (due to
standardized scales) while challenging DirectLiNGAM (due to heavy
measurement noise on the root variable, which triggers error
propagation). The true causal structure is a serial chain:
`x0 -> x1 -> x2 -> x3` (each coefficient 0.8).

## Usage

``` r
generate_lingam_paradox_data(n = 2000L, seed = 42L)
```

## Arguments

- n:

  number of samples (default: 2000)

- seed:

  random seed (default: 42)

## Value

list(data, true_adjacency)

- `data`: a data frame with 4 standardized variables (`x0`, `x1`, `x2`,
  `x3`); each column has a mean of 0 and a standard deviation of 1.

- `true_adjacency`: the 4x4 true adjacency matrix of the data-generating
  chain, following the `m[row = to, col = from]` convention and holding
  the structural coefficients (0.8) on the latent, pre-standardization
  scale.

## Details

This function intentionally injects strong measurement error into the
root (causal upstream) variable `x0`. This noise corrupts the
independence tests performed at the initial step of DirectLiNGAM,
frequently causing it to misidentify the root variable and leading to a
cascading failure (error propagation) throughout the causal ordering.

On the other hand, the output data is completely standardized using the
[`scale()`](https://rdrr.io/r/base/scale.html) function. This eliminates
any differences in scale among the variables, thereby neutralizing the
major weakness of ICA-LiNGAM (scale-dependence) and allowing it to
perform relatively better.

Because the data are standardized and the root carries measurement
error, the coefficients estimated by
[`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
will not exactly match the 0.8 values stored in `true_adjacency`.

## Examples

``` r
# Generate the dataset
paradox <- generate_lingam_paradox_data(n = 1000, seed = 123)

# Verify the dataset
head(paradox$data)
#>            x0        x1         x2         x3
#> 1 -0.98832138 -1.063113 -1.4745460 -1.8068974
#> 2  1.61196747  1.065200  0.1778246  1.0019652
#> 3  0.08506249 -0.909761 -1.3731009 -1.3846494
#> 4  1.75854498  1.845430  1.4870297  1.4747061
#> 5  0.47914188  2.009083  1.5757889  0.6862838
#> 6 -1.99217282 -1.410370 -0.8932110 -0.4404976
sapply(paradox$data, sd)
#> x0 x1 x2 x3 
#>  1  1  1  1 

# True data-generating structure
paradox$true_adjacency
#>     x0  x1  x2 x3
#> x0 0.0 0.0 0.0  0
#> x1 0.8 0.0 0.0  0
#> x2 0.0 0.8 0.0  0
#> x3 0.0 0.0 0.8  0
```
