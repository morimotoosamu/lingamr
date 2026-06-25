# plot QQ

plot QQ

## Usage

``` r
plot_residual_qq(X, lingam_result, ncol = 3, nrow = NULL)
```

## Arguments

- X:

  original data (matrix or data.frame)

- lingam_result:

  return value of lingam_direct()

- ncol:

  Number of columns.

- nrow:

  Number of rows.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object with QQ plots of residuals.

## Examples

``` r
# Load the sample data
LiNGAM_sample_1000 <- generate_lingam_sample_6()

# Run Direct LiNGAM
result <- lingam_direct(LiNGAM_sample_1000$data)

plot_residual_qq(LiNGAM_sample_1000$data, result)
```
