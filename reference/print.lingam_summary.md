# print method for lingam_summary

print method for lingam_summary

## Usage

``` r
# S3 method for class 'lingam_summary'
print(x, ...)
```

## Arguments

- x:

  A `lingam_summary` object

- ...:

  Additional arguments (unused)

## Value

The input object `x`, invisibly.

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()
model <- lingam_direct(LiNGAM_sample_1000$data, reg_method = "ols")
print(summary_lingam(LiNGAM_sample_1000$data, model))
#> === Direct LiNGAM Model Summary ===
#> Variables:    6
#> Observations: 1000
#> Edges:        15
#> Causal order: x3 -> x2 -> x0 -> x4 -> x5 -> x1
#> 
#> --- Assumption 1: Independence of residuals ---
#> Method:           spearman
#> Dependent pairs:  0 / 15  (p < 0.050)
#> Min p-value:      0.9187
#> => Residuals appear mutually independent (assumption supported).
#> 
#> --- Assumption 2: Non-Gaussianity of residuals ---
#> Method:           shapiro
#> Non-Gaussian:     6 / 6  (p <= 0.050)
#> => All residuals are non-Gaussian (assumption supported).
```
