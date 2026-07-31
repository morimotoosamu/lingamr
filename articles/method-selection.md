# Choosing a Causal Discovery Method

`lingamr` implements eleven causal discovery estimators. They all share
the same goal — recovering a directed causal structure from
observational data — but each relaxes a different assumption of the
basic LiNGAM model. This vignette helps you pick the right one.

The basic LiNGAM model assumes:

1.  **Linear** causal relationships,
2.  an **acyclic** causal graph (a DAG),
3.  **non-Gaussian**, mutually independent error terms,
4.  **no latent confounder** (every common cause is observed),
5.  **i.i.d.** observations (no time structure, no group structure),
6.  continuous variables and a complete (no `NA`) data matrix.

When all six hold, use
[`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md).
Every other estimator in the package exists to relax one (or two) of
these assumptions.

## Decision Guide

Work through the questions in order; the first “yes” points to the
method.

1.  **Are the observations a time series?**
    - Yes, and temporal dependence is autoregressive → **VAR-LiNGAM**
      ([`lingam_var()`](https://morimotoosamu.github.io/lingamr/reference/lingam_var.md))
    - Yes, and the VAR residuals stay autocorrelated even with more lags
      (moving-average disturbances) → **VARMA-LiNGAM**
      ([`lingam_varma()`](https://morimotoosamu.github.io/lingamr/reference/lingam_varma.md))
2.  **Are causal relationships plausibly nonlinear?**
    - Yes, and all major common causes are observed → **RESIT**
      ([`lingam_resit()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit.md))
    - Yes, and there may be unobserved variables → **CAM-UV**
      ([`lingam_camuv()`](https://morimotoosamu.github.io/lingamr/reference/lingam_camuv.md))
3.  **Might a latent confounder exist (linear case)?**
    - You want to know *which part of the causal order* is still
      trustworthy → **BottomUpParceLiNGAM**
      ([`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md))
    - You want to know *which specific pairs* are confounded → **RCD**
      ([`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md))
4.  **Does the data mix continuous and discrete variables?** → **LiM**
    ([`lingam_lim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_lim.md);
    binary by default, Poisson counts via `is_poisson = TRUE`)
5.  **Do you have several datasets that share a causal structure**
    (multiple sites, periods, cohorts)? → **MultiGroup Direct LiNGAM**
    ([`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md))
6.  **Does the data contain missing values (`NA`)?** →
    **[`bootstrap_with_imputation()`](https://morimotoosamu.github.io/lingamr/reference/bootstrap_with_imputation.md)**
    (bootstrap + multiple imputation)
7.  **Are there many variables (tens to hundreds), or even p \> n?** →
    **HighDimDirectLiNGAM**
    ([`lingam_high_dim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_high_dim.md))
8.  **None of the above** → **Direct LiNGAM**
    ([`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md))

If several complications apply at once (say, a nonlinear time series),
no single estimator handles both; you will have to prioritize the
violation that matters most for your data, or transform the data
(e.g. differencing a non-stationary series) so that fewer assumptions
are violated.

## Method Overview

| Method | Function | Handles | Output notes | Bootstrap |
|----|----|----|----|:--:|
| Direct LiNGAM | [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md) | The baseline: linear, acyclic, non-Gaussian, i.i.d. | Causal order + coefficient matrix | [`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md) |
| HighDimDirectLiNGAM | [`lingam_high_dim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_high_dim.md) | Many variables; p \> n | Same object as [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md) | — |
| VAR-LiNGAM | [`lingam_var()`](https://morimotoosamu.github.io/lingamr/reference/lingam_var.md) | Stationary time series (AR) | Instantaneous B0 + lagged matrices | [`lingam_var_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_var_bootstrap.md) |
| VARMA-LiNGAM | [`lingam_varma()`](https://morimotoosamu.github.io/lingamr/reference/lingam_varma.md) | Time series with MA errors | AR (psi) + MA (omega) matrices | [`lingam_varma_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_varma_bootstrap.md) |
| MultiGroup Direct LiNGAM | [`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md) | Multiple datasets, one shared order | Common order + per-group matrices | [`lingam_multi_group_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group_bootstrap.md) |
| BottomUpParceLiNGAM | [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md) | Latent confounders (linear) | Unresolved block; `NA` entries | [`lingam_parce_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce_bootstrap.md) |
| RCD | [`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md) | Latent confounders (linear) | Ancestor sets; confounded pairs `NA` | [`lingam_rcd_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd_bootstrap.md) |
| RESIT | [`lingam_resit()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit.md) | Nonlinear additive noise | 0/1 edges (no coefficients) | [`lingam_resit_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit_bootstrap.md) |
| CAM-UV | [`lingam_camuv()`](https://morimotoosamu.github.io/lingamr/reference/lingam_camuv.md) | Nonlinear + unobserved variables | Parents list; UCP/UBP pairs `NA` | — |
| LiM | [`lingam_lim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_lim.md) | Mixed continuous/discrete data | Coefficient matrix | — |
| Bootstrap with imputation | [`bootstrap_with_imputation()`](https://morimotoosamu.github.io/lingamr/reference/bootstrap_with_imputation.md) | Missing values (`NA`) | Aggregate via [`as_bootstrap_result()`](https://morimotoosamu.github.io/lingamr/reference/as_bootstrap_result.md) | (is itself a bootstrap) |

Two practical cost notes:

- **HSIC-based methods** (ParceLiNGAM, RCD, RESIT, CAM-UV, and
  `lingam_direct(measure = "kernel")`) build $`n \times n`$ Gram
  matrices per test; they are not recommended for $`n`$ in the thousands
  — subsample first.
- **Direct LiNGAM** costs $`O(p^3)`$ in the number of variables; for
  large $`p`$ switch to
  [`lingam_high_dim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_high_dim.md).

Worked examples for every method are on the package website:

- [Direct LiNGAM in
  depth](https://morimotoosamu.github.io/lingamr/articles/direct-lingam.html)
- [Bootstrap and
  diagnostics](https://morimotoosamu.github.io/lingamr/articles/bootstrap-diagnostics.html)
- [Time series: VAR and
  VARMA](https://morimotoosamu.github.io/lingamr/articles/time-series.html)
- [Latent confounders: ParceLiNGAM and
  RCD](https://morimotoosamu.github.io/lingamr/articles/latent-confounders.html)
- [Nonlinear methods: RESIT and
  CAM-UV](https://morimotoosamu.github.io/lingamr/articles/nonlinear.html)
- [Special data: LiM, MultiGroup, missing
  values](https://morimotoosamu.github.io/lingamr/articles/special-data.html)

## When LiNGAM Cannot Be Used

When the assumptions are not met, estimation either fails or
systematically recovers an incorrect structure.

| Assumption | When problems arise | Remedy / alternative |
|----|----|----|
| **Non-Gaussian errors** | When all errors follow a Gaussian distribution, the causal direction becomes unidentifiable | No LiNGAM variant can help; consider constraint-based methods (e.g. the PC algorithm in `pcalg`), which return an equivalence class instead of a unique direction |
| **Acyclic graph (DAG)** | When feedback loops (x -\> y -\> x) exist | Consider Cyclic LiNGAM (implemented in the Python `lingam` package) |
| **No latent common causes** | When unobserved common causes (hidden confounders) exist | [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md), [`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md) (linear), [`lingam_camuv()`](https://morimotoosamu.github.io/lingamr/reference/lingam_camuv.md) (nonlinear) |
| **Linear causal relationships** | When the relationships among variables are nonlinear | [`lingam_resit()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit.md), [`lingam_camuv()`](https://morimotoosamu.github.io/lingamr/reference/lingam_camuv.md) |
| **No measurement error (upstream variables)** | When heavy measurement error is present on variables near the root, the direction is systematically reversed | See the measurement error paradox in the [Direct LiNGAM article](https://morimotoosamu.github.io/lingamr/articles/direct-lingam.html) |
| **Independent and identically distributed (i.i.d.)** | Time-series data, or data pooled across heterogeneous sources | [`lingam_var()`](https://morimotoosamu.github.io/lingamr/reference/lingam_var.md) / [`lingam_varma()`](https://morimotoosamu.github.io/lingamr/reference/lingam_varma.md) (time series), [`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md) (grouped data) |
| **Sufficient sample size** | When $`n`$ is extremely small relative to the number of variables $`p`$ (rule of thumb: $`n < 10p`$), estimation tends to be unstable | Reduce the number of variables; sparsify with `reg_method = "adaptive_lasso"`; for $`p > n`$ use [`lingam_high_dim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_high_dim.md) |

## A Checklist to Verify in Advance

Before starting an actual analysis, we recommend confirming the
following.

1.  **Acyclicity of the graph** – Can feedback loops be ruled out from
    domain expertise?
2.  **Absence of latent variables** – Are the key observed variables all
    present? If not, prefer
    [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)
    /
    [`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)
    /
    [`lingam_camuv()`](https://morimotoosamu.github.io/lingamr/reference/lingam_camuv.md).
3.  **Non-Gaussianity of the errors** – Can be checked with
    [`test_residual_normality()`](https://morimotoosamu.github.io/lingamr/reference/test_residual_normality.md)
    (though this is a post-estimation diagnostic). As a quick check
    beforehand, visually inspect each variable’s histogram and skewness.
4.  **Presence of measurement error** – Is there measurement error on
    variables near the root? If so, interpret with care.
5.  **Sample size** – Aim for $`n \geq 10p`$. If it falls short, do not
    over-trust the results.

> **Summary:** LiNGAM is powerful when all assumptions – linear,
> acyclic, non-Gaussian, no latent variables, and i.i.d. – hold, and
> `lingamr` provides a dedicated estimator for the most common way each
> assumption breaks. Verifying the assumptions with domain knowledge and
> residual diagnostics before analysis is the first step toward reliable
> causal inference.
