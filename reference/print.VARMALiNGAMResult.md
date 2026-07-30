# Print method for VARMALiNGAMResult

Print method for VARMALiNGAMResult

## Usage

``` r
# S3 method for class 'VARMALiNGAMResult'
print(x, digits = 3, ...)
```

## Arguments

- x:

  VARMALiNGAMResult object

- digits:

  number of digits to display

- ...:

  additional arguments (unused)

## Value

The input object `x`, invisibly.

## Examples

``` r
sample <- generate_varmalingam_sample(n = 300, seed = 42)
model <- lingam_varma(sample$data,
  order = c(1, 1), criterion = NULL,
  reg_method = "ols", prune = FALSE
)
print(model)
#> VARMA-LiNGAM Result
#>   Variables : 3
#>   Order (p, q) : (1, 1)
#>   Causal order (instantaneous): x0 -> x1 -> x2
#> 
#> Instantaneous adjacency matrix B0 (row = to, col = from):
#>        x0    x1 x2
#> x0  0.000  0.00  0
#> x1  0.649  0.00  0
#> x2 -0.154 -0.42  0
#> 
#> Lagged adjacency matrix psi1 (row = to, col = from):
#>        x0     x1     x2
#> x0  0.412 -0.095  0.230
#> x1 -0.235  0.425 -0.084
#> x2  0.097  0.057  0.296
#> 
#> MA adjacency matrix omega1 (row = to, col = from):
#>        x0     x1     x2
#> x0  0.280  0.025 -0.131
#> x1 -0.022  0.124 -0.100
#> x2  0.058 -0.037  0.228
```
