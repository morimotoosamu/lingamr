# Summarize the goodness-of-fit of a Direct LiNGAM model at once

For a fitted Direct LiNGAM model, this verifies how well the two main
assumptions on which LiNGAM relies (mutual independence of residuals and
non-Gaussianity of residuals) hold, all at once, and displays the
results together. Internally it calls
[`get_error_independence_p_values()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values.md)
and
[`test_residual_normality()`](https://morimotoosamu.github.io/lingamr/reference/test_residual_normality.md).

## Usage

``` r
summary_lingam(
  X,
  lingam_result,
  independence_method = "spearman",
  normality_method = "shapiro",
  alpha = 0.05
)
```

## Arguments

- X:

  The original data (matrix or data.frame), the one used to estimate
  `lingam_result`.

- lingam_result:

  The return value of
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
  (a `LingamResult` object)

- independence_method:

  The type of correlation coefficient used in the residual independence
  test ("spearman", "pearson", "kendall"). Passed to
  [`get_error_independence_p_values()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values.md).

- normality_method:

  The method for the residual normality test ("shapiro", "ks", "ad",
  "lillie", "jb"). Passed to
  [`test_residual_normality()`](https://morimotoosamu.github.io/lingamr/reference/test_residual_normality.md).

- alpha:

  Significance level (default: 0.05)

## Value

A list of class `lingam_summary`, containing the following elements:

- `n_variables`, `n_samples`: Number of variables / number of
  observations

- `causal_order`: Causal order (variable-name labels)

- `n_edges`: Number of nonzero elements in the adjacency matrix (number
  of estimated edges)

- `independence_p_values`: Matrix of p-values from the independence test
  between residuals

- `n_dependent_pairs`, `n_pairs`: Number of pairs with p \< alpha /
  total number of pairs

- `min_independence_p`: Minimum p-value of the independence test

- `normality`: Result of the normality test (a `lingam_normality_test`
  object)

- `n_non_gaussian`: Number of variables judged to be non-Gaussian

- `alpha`, `independence_method`, `normality_method`: The settings used

## Details

Gaussian-likelihood-based criteria such as BIC/AIC are not included
because they are theoretically inconsistent with LiNGAM's assumption
that "the errors are non-Gaussian". Instead, the verification results of
the assumptions themselves are summarized.

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

model <- lingam_direct(LiNGAM_sample_1000$data, reg_method = "ols")

summary_lingam(LiNGAM_sample_1000$data, model)
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
