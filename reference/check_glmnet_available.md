# Check whether glmnet is available

If it is not available, raise an error indicating which regression
method required it.

## Usage

``` r
check_glmnet_available(method)
```

## Arguments

- method:

  name of the regression method that requires glmnet (for the error
  message)
