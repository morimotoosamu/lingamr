# Resolve the regressor argument to a fitting function and its label

The label is normalized (`match.arg` accepts partial matches such as
`"g"`), so the `regressor` field of the result always reports the
canonical name of what actually ran.

## Usage

``` r
resit_make_regressor(regressor)
```

## Arguments

- regressor:

  either a string (currently only `"gam"`) or a function
  `function(X, y)` returning fitted values

## Value

`list(fn = <function (X, y) -> fitted values>, label = <string>)`
