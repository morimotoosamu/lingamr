## Resubmission

This is a resubmission. In this version I have:

* Put software names in single quotes in Title and Description
  (e.g., 'LiNGAM', 'Python', 'VAR-LiNGAM').
* Added `\value` tags to all exported function documentation
  (plot_residual_qq, print.BootstrapResult, print.causal_order_stability,
  print.lingam_normality_test, print.lingam_summary, print.LingamResult).

## R CMD check results

0 errors | 0 warnings | 1 note

The NOTE is from "checking R code for possible problems":

* autoplot.LingamResult: no visible binding for global variable 'lx' / 'ly'.
  These are column names used inside ggplot2::aes() and are standard
  non-standard evaluation usage.

## Test environments

* Local: Windows 11, R 4.6.0

## Downstream dependencies

There are currently no downstream dependencies for this package.
