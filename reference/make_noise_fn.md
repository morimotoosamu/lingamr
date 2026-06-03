# Create noise generation function

Internal helper to create a noise function for the specified
distribution.

## Usage

``` r
make_noise_fn(noise_dist)
```

## Arguments

- noise_dist:

  distribution name "uniform" : Uniform(0, 1) - non-Gaussian (LiNGAM
  works well) "gaussian" : Normal(0, 1) - LiNGAM may fail "lognormal" :
  Log-normal(0, 1) - skewed, non-Gaussian "exponential" :
  Exponential(1) - skewed, non-Gaussian "t3" : t-distribution (df=3) -
  heavy tails

## Value

function(n) that generates n random numbers
