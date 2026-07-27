# Emit the shared bootstrap completion message

Emit the shared bootstrap completion message

## Usage

``` r
bootstrap_completion_message(t_start, n_success, n_sampling)
```

## Arguments

- t_start:

  Result of [`proc.time()`](https://rdrr.io/r/base/proc.time.html)
  captured before the iterations began.

- n_success:

  Number of iterations that succeeded.

- n_sampling:

  Number of iterations that were requested.

## Value

`NULL`, invisibly. Called for the
[`message()`](https://rdrr.io/r/base/message.html) side effect.
