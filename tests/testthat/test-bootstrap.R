test_that("lingam_direct_bootstrap validates inputs before running", {
  dat <- generate_lingam_sample_6(n = 100, seed = 1)

  expect_error(lingam_direct_bootstrap(dat$data, n_sampling = 5L, measure = "bad"))
  expect_error(lingam_direct_bootstrap(dat$data, n_sampling = 5L, reg_method = "bad"))
  expect_error(lingam_direct_bootstrap(dat$data, n_sampling = 5L, lambda = "bad"))
  expect_error(lingam_direct_bootstrap(dat$data, n_sampling = 5L, init_method = "bad"))
  expect_error(
    lingam_direct_bootstrap(dat$data, n_sampling = "abc"),
    "n_sampling must be a positive integer"
  )
  expect_error(
    lingam_direct_bootstrap(dat$data, n_sampling = 0L),
    "n_sampling must be a positive integer"
  )
  expect_error(
    lingam_direct_bootstrap(dat$data, n_sampling = 5L, compute_total_effects = "yes"),
    "compute_total_effects must be"
  )
})

test_that("compute_total_effects = FALSE skips total effects but keeps edge stability intact", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)

  bs <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, reg_method = "ols",
                                 seed = 1, compute_total_effects = FALSE)

  expect_null(bs$total_effects)
  expect_true(is.array(bs$adjacency_matrices))
  expect_equal(dim(bs$adjacency_matrices)[1], 10L)
  expect_equal(dim(get_probabilities(bs)), c(6L, 6L))
  expect_error(
    get_total_causal_effects(bs),
    "compute_total_effects = FALSE"
  )
})

test_that("lingam_direct_bootstrap passes init_method through", {
  skip_if_not_installed("glmnet")
  dat <- generate_lingam_sample_6(n = 100, seed = 1)

  bs <- lingam_direct_bootstrap(dat$data,
    n_sampling = 2L, seed = 42L, verbose = FALSE,
    reg_method = "adaptive_lasso", init_method = "ridge"
  )
  expect_s3_class(bs, "BootstrapResult")
})

test_that("lingam_direct_bootstrap returns BootstrapResult", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, reg_method = "ols", seed = 42L)

  expect_s3_class(bs, "BootstrapResult")
  expect_true(is.array(bs$adjacency_matrices))
  expect_equal(dim(bs$adjacency_matrices)[1], 10L)
  expect_equal(dim(bs$adjacency_matrices)[2], 6L)
  expect_equal(dim(bs$adjacency_matrices)[3], 6L)
})

test_that("lingam_direct_bootstrap is reproducible with same seed", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)

  bs1 <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, reg_method = "ols", seed = 99L)
  bs2 <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, reg_method = "ols", seed = 99L)

  expect_equal(bs1$adjacency_matrices, bs2$adjacency_matrices)
})

test_that("lingam_direct_bootstrap with different seeds gives different results", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)

  bs1 <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, reg_method = "ols", seed = 1L)
  bs2 <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, reg_method = "ols", seed = 2L)

  expect_false(identical(bs1$adjacency_matrices, bs2$adjacency_matrices))
})

test_that("get_causal_direction_counts returns data.frame with expected columns", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 15L, reg_method = "ols", seed = 42L)
  dc  <- get_causal_direction_counts(bs)

  expect_s3_class(dc, "data.frame")
  expect_true(all(c("from", "to", "count", "proportion",
                    "mean_effect", "median_effect", "sd_effect",
                    "ci_lower", "ci_upper") %in% names(dc)))
  # from comes before to (column order)
  expect_lt(which(names(dc) == "from"), which(names(dc) == "to"))
  # proportion is in [0, 1]
  expect_true(all(dc$proportion >= 0 & dc$proportion <= 1))
})

test_that("get_causal_direction_counts with split_by_causal_effect_sign adds sign column", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 15L, reg_method = "ols", seed = 42L)
  dc  <- get_causal_direction_counts(bs, split_by_causal_effect_sign = TRUE)

  expect_true("sign" %in% names(dc))
  expect_true(all(dc$sign %in% c(-1L, 1L)))
})

test_that("get_adjacency_matrix_summary returns correctly shaped matrix", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 15L, reg_method = "ols", seed = 42L)
  B   <- get_adjacency_matrix_summary(bs)

  expect_true(is.matrix(B))
  expect_equal(dim(B), c(6L, 6L))
})

# =============================================================================
# Reproducibility tests for parallel execution (#11)
# =============================================================================

test_that("parallel bootstrap is reproducible with same seed and same n_cores", {
  skip_on_cran()
  skip_if_not(parallel::detectCores() >= 2L, "requires >= 2 cores")

  dat <- generate_lingam_sample_6(n = 300, seed = 1)

  bs1 <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, reg_method = "ols", seed = 77L,
                                  parallel = TRUE, n_cores = 2L)
  bs2 <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, reg_method = "ols", seed = 77L,
                                  parallel = TRUE, n_cores = 2L)

  expect_equal(bs1$adjacency_matrices, bs2$adjacency_matrices)
})

test_that("parallel bootstrap with different seeds gives different results", {
  skip_on_cran()
  skip_if_not(parallel::detectCores() >= 2L, "requires >= 2 cores")

  dat <- generate_lingam_sample_6(n = 300, seed = 1)

  bs1 <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, reg_method = "ols", seed = 10L,
                                  parallel = TRUE, n_cores = 2L)
  bs2 <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, reg_method = "ols", seed = 20L,
                                  parallel = TRUE, n_cores = 2L)

  expect_false(identical(bs1$adjacency_matrices, bs2$adjacency_matrices))
})

test_that("parallel and sequential results differ (L'Ecuyer vs set.seed)", {
  # As documented, parallel and sequential runs do not produce identical
  # numbers even with the same seed.
  # This test records that specification as a regression test.
  skip_on_cran()
  skip_if_not(parallel::detectCores() >= 2L, "requires >= 2 cores")

  dat <- generate_lingam_sample_6(n = 300, seed = 1)

  bs_seq <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, reg_method = "ols", seed = 42L,
                                    parallel = FALSE)
  bs_par <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, reg_method = "ols", seed = 42L,
                                    parallel = TRUE, n_cores = 2L)

  expect_false(identical(bs_seq$adjacency_matrices, bs_par$adjacency_matrices))
})

test_that("parallel bootstrap returns same structure as sequential", {
  skip_on_cran()
  skip_if_not(parallel::detectCores() >= 2L, "requires >= 2 cores")

  dat <- generate_lingam_sample_6(n = 300, seed = 1)

  bs_seq <- lingam_direct_bootstrap(dat$data, n_sampling = 8L, reg_method = "ols", seed = 1L,
                                    parallel = FALSE)
  bs_par <- lingam_direct_bootstrap(dat$data, n_sampling = 8L, reg_method = "ols", seed = 1L,
                                    parallel = TRUE, n_cores = 2L)

  # Even if the numbers differ, the structure (dimensions/class) must be the same
  expect_s3_class(bs_par, "BootstrapResult")
  expect_equal(dim(bs_par$adjacency_matrices), dim(bs_seq$adjacency_matrices))
})

# =============================================================================
# Value tests for the BootstrapResult query helpers (previously untested;
# their VAR-side counterparts already have equivalent value tests)
# =============================================================================

test_that("get_probabilities has the right shape and detects the true edges", {
  dat <- generate_lingam_sample_6(n = 800, seed = 42)
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 20L, reg_method = "ols", seed = 1)

  p <- get_probabilities(bs, min_causal_effect = 0.3)
  expect_equal(dim(p), c(6L, 6L))
  expect_true(all(p >= 0 & p <= 1))
  # x3 -> x0 (true coefficient 3.0) should be detected in nearly every resample
  expect_gt(p[1, 4], 0.8)
  # x1 has no children in the true DAG, so any edge out of it (column 2) should
  # only ever show up as small OLS estimation noise below a reasonable threshold
  expect_true(all(p[, 2] == 0))
})

test_that("get_paths finds the indirect path x3 -> x0 -> x1", {
  dat <- generate_lingam_sample_6(n = 800, seed = 42)
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 20L, reg_method = "ols", seed = 1)

  paths <- get_paths(bs, from_index = 4, to_index = 2)
  expect_s3_class(paths, "data.frame")
  expect_named(paths, c("path", "effect", "probability"))
  expect_gt(nrow(paths), 0)
  expect_true(all(paths$probability > 0 & paths$probability <= 1))
  has_chain <- any(vapply(paths$path, function(p) identical(p, c(4L, 1L, 2L)), logical(1)))
  expect_true(has_chain)
})

test_that("get_total_causal_effects reports the known strong edges", {
  dat <- generate_lingam_sample_6(n = 800, seed = 42)
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 20L, reg_method = "ols", seed = 1)

  te <- get_total_causal_effects(bs, min_causal_effect = 0.3)
  expect_s3_class(te, "data.frame")
  expect_named(te, c("from", "to", "effect", "probability"))
  expect_true(all(te$probability > 0 & te$probability <= 1))
  # x3 -> x0 direct effect (true 3.0)
  row <- te[te$from == 4 & te$to == 1, ]
  expect_equal(nrow(row), 1L)
  expect_equal(row$effect, 3.0, tolerance = 0.5)
})

test_that("plot_bootstrap_probabilities returns a grViz object", {
  skip_if_not_installed("DiagrammeR")
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, reg_method = "ols", seed = 1)

  g <- plot_bootstrap_probabilities(bs)
  expect_s3_class(g, "grViz")
})

test_that("get_directed_acyclic_graph_counts returns dag/count lists that sum to n_sampling", {
  dat <- generate_lingam_sample_6(n = 500, seed = 42)
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 20L, reg_method = "ols", seed = 1)

  res <- get_directed_acyclic_graph_counts(bs)
  expect_named(res, c("dag", "count"))
  expect_true(is.list(res$dag))
  expect_equal(length(res$dag), length(res$count))
  expect_equal(sum(res$count), 20L)
  # counts are sorted in descending order
  expect_true(all(diff(res$count) <= 0))
  expect_true(all(vapply(res$dag, is.data.frame, logical(1))))
})
