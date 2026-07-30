# Package index

## Direct LiNGAM

Core causal discovery on i.i.d. data

- [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
  : Direct LiNGAM
- [`lingam_high_dim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_high_dim.md)
  : High-Dimensional Direct LiNGAM
- [`make_prior_knowledge()`](https://morimotoosamu.github.io/lingamr/reference/make_prior_knowledge.md)
  : Create a prior knowledge matrix
- [`summary_lingam()`](https://morimotoosamu.github.io/lingamr/reference/summary_lingam.md)
  : Summarize the goodness-of-fit of a Direct LiNGAM model at once
- [`estimate_total_effect()`](https://morimotoosamu.github.io/lingamr/reference/estimate_total_effect.md)
  : Estimate the total causal effect between two specified variables
- [`estimate_all_total_effects()`](https://morimotoosamu.github.io/lingamr/reference/estimate_all_total_effects.md)
  : Estimate the total causal effects between all variables at once
- [`get_error_independence_p_values()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values.md)
  : Compute p-values for the independence test of the errors
- [`get_causal_order_stability()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_order_stability.md)
  : Evaluate the stability of the causal order from bootstrap

## Bootstrap stability analysis

- [`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md)
  : Bootstrap for Direct LiNGAM
- [`get_probabilities()`](https://morimotoosamu.github.io/lingamr/reference/get_probabilities.md)
  : Get bootstrap probabilities
- [`get_causal_direction_counts()`](https://morimotoosamu.github.io/lingamr/reference/get_causal_direction_counts.md)
  : Get counts, proportions, and causal effects of causal directions
- [`get_directed_acyclic_graph_counts()`](https://morimotoosamu.github.io/lingamr/reference/get_directed_acyclic_graph_counts.md)
  : Get DAG counts
- [`get_adjacency_matrix_summary()`](https://morimotoosamu.github.io/lingamr/reference/get_adjacency_matrix_summary.md)
  : Create an adjacency matrix of representative causal-effect values
  from bootstrap results
- [`get_total_causal_effects()`](https://morimotoosamu.github.io/lingamr/reference/get_total_causal_effects.md)
  : Get a list of total causal effects
- [`get_paths()`](https://morimotoosamu.github.io/lingamr/reference/get_paths.md)
  : Get all paths between two specified variables and their bootstrap
  probabilities
- [`plot_bootstrap_probabilities()`](https://morimotoosamu.github.io/lingamr/reference/plot_bootstrap_probabilities.md)
  : Draw bootstrap probabilities with DiagrammeR

## VAR-LiNGAM (time series)

- [`lingam_var()`](https://morimotoosamu.github.io/lingamr/reference/lingam_var.md)
  : VAR-LiNGAM for time series causal discovery
- [`lingam_var_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_var_bootstrap.md)
  : Bootstrap for VAR-LiNGAM
- [`get_var_probabilities()`](https://morimotoosamu.github.io/lingamr/reference/get_var_probabilities.md)
  : Bootstrap probabilities for a VAR-LiNGAM model
- [`get_var_paths()`](https://morimotoosamu.github.io/lingamr/reference/get_var_paths.md)
  : Enumerate bootstrap paths between two variables in a VAR-LiNGAM
  model
- [`estimate_var_total_effect()`](https://morimotoosamu.github.io/lingamr/reference/estimate_var_total_effect.md)
  : Estimate a total causal effect in a VAR-LiNGAM model
- [`check_var_stationarity()`](https://morimotoosamu.github.io/lingamr/reference/check_var_stationarity.md)
  : Check the stationarity of a fitted VAR-LiNGAM model
- [`test_varlingam_residual_normality()`](https://morimotoosamu.github.io/lingamr/reference/test_varlingam_residual_normality.md)
  : Test the non-Gaussianity of VAR-LiNGAM residuals
- [`test_varlingam_residual_normality_all()`](https://morimotoosamu.github.io/lingamr/reference/test_varlingam_residual_normality_all.md)
  : Run several normality tests on VAR-LiNGAM residuals at once
- [`plot_varlingam_residual_qq()`](https://morimotoosamu.github.io/lingamr/reference/plot_varlingam_residual_qq.md)
  : Q-Q plots of VAR-LiNGAM residuals

## MultiGroup Direct LiNGAM

Joint estimation across multiple datasets sharing a causal order

- [`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md)
  : Multi-Group Direct LiNGAM
- [`lingam_multi_group_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group_bootstrap.md)
  : Bootstrap for Multi-Group Direct LiNGAM
- [`get_group_result()`](https://morimotoosamu.github.io/lingamr/reference/get_group_result.md)
  : Extract a single group's result from a MultiGroupLingamResult

## BottomUpParceLiNGAM

Causal discovery robust against latent confounders

- [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)
  : Bottom-Up ParceLiNGAM
- [`lingam_parce_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce_bootstrap.md)
  : Bootstrap for Bottom-Up ParceLiNGAM
- [`estimate_total_effect_parce()`](https://morimotoosamu.github.io/lingamr/reference/estimate_total_effect_parce.md)
  : Estimate the total causal effect between two variables (ParceLiNGAM)
- [`get_error_independence_p_values_parce()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values_parce.md)
  : Compute p-values for the independence of ParceLiNGAM residuals
  (HSIC-based)

## RCD (Repetitive Causal Discovery)

Causal discovery with latent confounder detection

- [`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)
  : RCD (Repetitive Causal Discovery)
- [`lingam_rcd_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd_bootstrap.md)
  : Bootstrap for RCD
- [`estimate_total_effect_rcd()`](https://morimotoosamu.github.io/lingamr/reference/estimate_total_effect_rcd.md)
  : Estimate the total causal effect between two variables (RCD)
- [`get_error_independence_p_values_rcd()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values_rcd.md)
  : Compute p-values for the independence of RCD residuals (HSIC-based)

## RESIT (nonlinear additive noise models)

- [`lingam_resit()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit.md)
  : RESIT causal discovery for nonlinear additive noise models
- [`lingam_resit_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit_bootstrap.md)
  : Bootstrap for RESIT

## CAM-UV (nonlinear models with unobserved variables)

- [`lingam_camuv()`](https://morimotoosamu.github.io/lingamr/reference/lingam_camuv.md)
  : CAM-UV (Causal Additive Models with Unobserved Variables)

## LiM (mixed continuous and discrete data)

- [`lingam_lim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_lim.md)
  : LiM: LiNGAM for Mixed Data

## Missing data

- [`bootstrap_with_imputation()`](https://morimotoosamu.github.io/lingamr/reference/bootstrap_with_imputation.md)
  : Bootstrap with Multiple Imputation for Direct LiNGAM
- [`as_bootstrap_result()`](https://morimotoosamu.github.io/lingamr/reference/as_bootstrap_result.md)
  : Collapse an ImputationBootstrapResult into a BootstrapResult

## Model fit and diagnostics

- [`evaluate_model_fit()`](https://morimotoosamu.github.io/lingamr/reference/evaluate_model_fit.md)
  : Evaluate model fit of an estimated causal graph
- [`test_residual_normality()`](https://morimotoosamu.github.io/lingamr/reference/test_residual_normality.md)
  : Test normality of residuals from Direct LiNGAM
- [`plot_residual_qq()`](https://morimotoosamu.github.io/lingamr/reference/plot_residual_qq.md)
  : plot QQ
- [`plot_adjacency()`](https://morimotoosamu.github.io/lingamr/reference/plot_adjacency.md)
  : Plot a causal graph from an adjacency matrix with DiagrammeR

## Sample data generation

- [`generate_camuv_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_camuv_sample.md)
  : Generate sample data with unobserved variables (for CAM-UV)
- [`generate_lim_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_lim_sample.md)
  : Generate sample data for LiM (3 mixed variables)
- [`generate_lingam_hard_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_lingam_hard_sample.md)
  : Generate a challenging sample data for Direct LiNGAM
- [`generate_lingam_large_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_lingam_large_sample.md)
  : Generate large-scale sample data to benchmark Direct LiNGAM
  scalability
- [`generate_lingam_paradox_data()`](https://morimotoosamu.github.io/lingamr/reference/generate_lingam_paradox_data.md)
  : Generate Paradoxical Data Where DirectLiNGAM Struggles
- [`generate_lingam_sample_10()`](https://morimotoosamu.github.io/lingamr/reference/generate_lingam_sample_10.md)
  : Generate 10-variable sample data for Direct LiNGAM
- [`generate_lingam_sample_6()`](https://morimotoosamu.github.io/lingamr/reference/generate_lingam_sample_6.md)
  : Generate sample data for Direct LiNGAM (6 variables)
- [`generate_multi_group_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_multi_group_sample.md)
  : Generate sample data for Multi-Group Direct LiNGAM (2 groups, 6
  variables)
- [`generate_parce_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_parce_sample.md)
  : Generate sample data with a latent confounder (for
  BottomUpParceLiNGAM)
- [`generate_rcd_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_rcd_sample.md)
  : Generate sample data with a latent confounder (for RCD)
- [`generate_resit_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_resit_sample.md)
  : Generate sample data from a nonlinear additive noise model (for
  RESIT)
- [`generate_varlingam_sample()`](https://morimotoosamu.github.io/lingamr/reference/generate_varlingam_sample.md)
  : Generate sample data from a VAR-LiNGAM model

## Tidiers and autoplot

- [`tidy(`*`<BootstrapResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/tidy.BootstrapResult.md)
  : Convert a BootstrapResult to a tidy data.frame
- [`tidy(`*`<CAMUVResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/tidy.CAMUVResult.md)
  : Convert a CAMUVResult to a tidy data.frame
- [`tidy(`*`<ImputationBootstrapResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/tidy.ImputationBootstrapResult.md)
  : Convert an ImputationBootstrapResult to a tidy data.frame
- [`tidy(`*`<LiMResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/tidy.LiMResult.md)
  : Convert a LiMResult to a tidy data.frame
- [`tidy(`*`<LingamResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/tidy.LingamResult.md)
  : Convert a LingamResult to a tidy data.frame
- [`tidy(`*`<MultiGroupBootstrapResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/tidy.MultiGroupBootstrapResult.md)
  : Convert a MultiGroupBootstrapResult to a tidy data.frame
- [`tidy(`*`<MultiGroupLingamResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/tidy.MultiGroupLingamResult.md)
  : Convert a MultiGroupLingamResult to a tidy data.frame
- [`tidy(`*`<ParceLingamResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/tidy.ParceLingamResult.md)
  : Convert a ParceLingamResult to a tidy data.frame
- [`tidy(`*`<RCDResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/tidy.RCDResult.md)
  : Convert an RCDResult to a tidy data.frame
- [`tidy(`*`<ResitResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/tidy.ResitResult.md)
  : Convert a ResitResult to a tidy data.frame
- [`glance(`*`<CAMUVResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/glance.CAMUVResult.md)
  : Get a one-row summary of a CAMUVResult
- [`glance(`*`<LiMResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/glance.LiMResult.md)
  : Get a one-row summary of a LiMResult
- [`glance(`*`<LingamResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/glance.LingamResult.md)
  : Get a one-row summary of a LingamResult
- [`glance(`*`<MultiGroupLingamResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/glance.MultiGroupLingamResult.md)
  : Get a one-row summary of a MultiGroupLingamResult
- [`glance(`*`<ParceLingamResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/glance.ParceLingamResult.md)
  : Get a one-row summary of a ParceLingamResult
- [`glance(`*`<RCDResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/glance.RCDResult.md)
  : Get a one-row summary of an RCDResult
- [`glance(`*`<ResitResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/glance.ResitResult.md)
  : Get a one-row summary of a ResitResult
- [`autoplot(`*`<CAMUVResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/autoplot.CAMUVResult.md)
  : Plot the causal graph of a CAMUVResult with ggplot2
- [`autoplot(`*`<LiMResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/autoplot.LiMResult.md)
  : Plot the causal graph of a LiMResult with ggplot2
- [`autoplot(`*`<LingamResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/autoplot.LingamResult.md)
  : Plot the causal graph of a LingamResult with ggplot2
- [`autoplot(`*`<MultiGroupLingamResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/autoplot.MultiGroupLingamResult.md)
  : Plot one group of a MultiGroupLingamResult with ggplot2
- [`autoplot(`*`<ParceLingamResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/autoplot.ParceLingamResult.md)
  : Plot the causal graph of a ParceLingamResult with ggplot2
- [`autoplot(`*`<RCDResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/autoplot.RCDResult.md)
  : Plot the causal graph of an RCDResult with ggplot2
- [`autoplot(`*`<ResitResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/autoplot.ResitResult.md)
  : Plot the causal graph of a ResitResult with ggplot2

## Print methods

- [`print(`*`<BootstrapResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/print.BootstrapResult.md)
  : Display the contents of a BootstrapResult
- [`print(`*`<CAMUVResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/print.CAMUVResult.md)
  : Print method for CAMUVResult
- [`print(`*`<ImputationBootstrapResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/print.ImputationBootstrapResult.md)
  : Print method for ImputationBootstrapResult
- [`print(`*`<LiMResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/print.LiMResult.md)
  : Print method for LiMResult
- [`print(`*`<LingamResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/print.LingamResult.md)
  : Print method for LingamResult
- [`print(`*`<MultiGroupBootstrapResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/print.MultiGroupBootstrapResult.md)
  : Print method for MultiGroupBootstrapResult
- [`print(`*`<MultiGroupLingamResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/print.MultiGroupLingamResult.md)
  : Print method for MultiGroupLingamResult
- [`print(`*`<ParceLingamResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/print.ParceLingamResult.md)
  : Print method for ParceLingamResult
- [`print(`*`<RCDResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/print.RCDResult.md)
  : Print method for RCDResult
- [`print(`*`<ResitResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/print.ResitResult.md)
  : Print method for ResitResult
- [`print(`*`<VARBootstrapResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/print.VARBootstrapResult.md)
  : Print a VARBootstrapResult
- [`print(`*`<VARLiNGAMResult>`*`)`](https://morimotoosamu.github.io/lingamr/reference/print.VARLiNGAMResult.md)
  : Print method for VARLiNGAMResult
- [`print(`*`<causal_order_stability>`*`)`](https://morimotoosamu.github.io/lingamr/reference/print.causal_order_stability.md)
  : print method for causal_order_stability
- [`print(`*`<lingam_normality_test>`*`)`](https://morimotoosamu.github.io/lingamr/reference/print.lingam_normality_test.md)
  : Print method for lingam_normality_test
- [`print(`*`<lingam_summary>`*`)`](https://morimotoosamu.github.io/lingamr/reference/print.lingam_summary.md)
  : print method for lingam_summary
- [`print(`*`<var_stationarity>`*`)`](https://morimotoosamu.github.io/lingamr/reference/print.var_stationarity.md)
  : Print method for var_stationarity
