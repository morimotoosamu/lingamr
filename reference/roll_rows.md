# Roll matrix rows (numpy np.roll equivalent, axis = 0)

Shifts the rows of `M` downward by `shift`, wrapping the last `shift`
rows around to the top. Used to build the lagged design for total-effect
regression. The wrap-around contaminates the first `shift` rows,
matching the Python reference (the effect is negligible for long
series).

## Usage

``` r
roll_rows(M, shift)
```

## Arguments

- M:

  numeric matrix

- shift:

  non-negative integer number of rows to shift down

## Value

matrix with rolled rows
