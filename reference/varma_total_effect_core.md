# Total-effect regression core for VARMA-LiNGAM

Regresses the destination variable on the source variable plus the
source's parents (a back-door adjustment) over a joined design of lagged
X and lagged LiNGAM residuals, and returns the source's coefficient.
Shared by
[`estimate_varma_total_effect()`](https://morimotoosamu.github.io/lingamr/reference/estimate_varma_total_effect.md)
and the VARMA-LiNGAM bootstrap.

## Usage

``` r
varma_total_effect_core(
  X,
  ee_full,
  am_joined,
  order,
  from_index,
  to_index,
  from_lag,
  X_joined = NULL
)
```

## Arguments

- X:

  numeric matrix (n x m), rows ordered in time

- ee_full:

  LiNGAM residuals `e_t = (I - B0) n_t`, full length (n x m)

- am_joined:

  joined causal matrix `cbind(psi_0..psi_p, omega_1..omega_q)` (m x
  m(1 + p + q))

- order:

  VARMA order c(p, q)

- from_index:

  source variable (1-based integer)

- to_index:

  destination variable (1-based integer)

- from_lag:

  lag of the source variable (non-negative integer)

- X_joined:

  precomputed
  [`varma_joined_design()`](https://morimotoosamu.github.io/lingamr/reference/varma_joined_design.md)
  result, or NULL to build it here. The design depends only on (X,
  ee_full, order, from_lag), so the bootstrap builds it once per lag
  instead of once per variable pair.

## Value

the estimated total effect (scalar)

## Details

The X region of the design holds lags `0..(p + from_lag)` and the
residual region holds lags `1..(q + from_lag)`, so parents shifted into
the `from_lag` block always reference filled columns. (The Python
reference allocates the extra `from_lag` blocks but never fills them, so
its shifted predictors can point at zero or misaligned columns when
`from_lag > 0`.)
