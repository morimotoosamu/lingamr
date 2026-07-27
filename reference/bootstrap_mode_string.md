# Human-readable execution-mode label for bootstrap progress messages

Human-readable execution-mode label for bootstrap progress messages

## Usage

``` r
bootstrap_mode_string(parallel, n_cores)
```

## Arguments

- parallel:

  Whether the run is parallel (logical, after core resolution).

- n_cores:

  Resolved core count (used only when `parallel` is `TRUE`).

## Value

`"parallel, N cores"` or `"sequential"`.
