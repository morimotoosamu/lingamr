# Generate 10-variable sample data for Direct LiNGAM

Generates a sample dataset with a known causal structure. The true
causal structure is: x3 -\> x0 (coef = 3.0) x3 -\> x2 (coef = 6.0) x3
-\> x9 (coef = 7.0) x0 -\> x1 (coef = 3.0) x0 -\> x5 (coef = 4.0) x0 -\>
x4 (coef = 8.0) x0 -\> x7 (coef = 3.0) x2 -\> x1 (coef = 2.0) x2 -\> x4
(coef = -1.0) x2 -\> x8 (coef = 0.5) x1 -\> x6 (coef = 2.0) x5 -\> x8
(coef = 2.0) x4 -\> x7 (coef = 1.5) x6 -\> x9 (coef = 1.0)

## Usage

``` r
generate_lingam_sample_10(n = 1000L, seed = 42L, noise_dist = "uniform")
```

## Arguments

- n:

  number of samples (default: 1000)

- seed:

  random seed (default: 42)

- noise_dist:

  error term distribution "uniform" : Uniform(0, 1) - default,
  non-Gaussian (LiNGAM works well) "gaussian" : Normal(0, 1) - LiNGAM
  may fail "lognormal" : Log-normal(0, 1) - skewed, non-Gaussian
  "exponential" : Exponential(1) - skewed, non-Gaussian "t3" :
  t-distribution (df=3) - heavy tails

## Value

list(data, true_adjacency)

## Examples

``` r
# Non-Gaussian (LiNGAM works well)
X_nongauss <- generate_lingam_sample_10(noise_dist = "uniform")
result <- lingam_direct(X_nongauss$data)
result$causal_order
#>  [1]  4  3  1  5  6  8  2  9  7 10

# Gaussian (LiNGAM may fail)
X_gauss <- generate_lingam_sample_10(noise_dist = "gaussian")
result <- lingam_direct(X_gauss$data)
result$causal_order
#>  [1] 10  2  7  3  9  4  6  5  8  1
```
