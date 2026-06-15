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
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, seed = 42L)

  expect_s3_class(bs, "BootstrapResult")
  expect_true(is.array(bs$adjacency_matrices))
  expect_equal(dim(bs$adjacency_matrices)[1], 10L)
  expect_equal(dim(bs$adjacency_matrices)[2], 6L)
  expect_equal(dim(bs$adjacency_matrices)[3], 6L)
})

test_that("lingam_direct_bootstrap is reproducible with same seed", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)

  bs1 <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, seed = 99L)
  bs2 <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, seed = 99L)

  expect_equal(bs1$adjacency_matrices, bs2$adjacency_matrices)
})

test_that("lingam_direct_bootstrap with different seeds gives different results", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)

  bs1 <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, seed = 1L)
  bs2 <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, seed = 2L)

  expect_false(identical(bs1$adjacency_matrices, bs2$adjacency_matrices))
})

test_that("get_causal_direction_counts returns data.frame with expected columns", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 15L, seed = 42L)
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
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 15L, seed = 42L)
  dc  <- get_causal_direction_counts(bs, split_by_causal_effect_sign = TRUE)

  expect_true("sign" %in% names(dc))
  expect_true(all(dc$sign %in% c(-1L, 1L)))
})

test_that("get_adjacency_matrix_summary returns correctly shaped matrix", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 15L, seed = 42L)
  B   <- get_adjacency_matrix_summary(bs)

  expect_true(is.matrix(B))
  expect_equal(dim(B), c(6L, 6L))
})

# =============================================================================
# Reproducibility tests for parallel execution (#11)
# =============================================================================

test_that("parallel bootstrap is reproducible with same seed and same n_cores", {
  skip_if_not(parallel::detectCores() >= 2L, "requires >= 2 cores")

  dat <- generate_lingam_sample_6(n = 300, seed = 1)

  bs1 <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, seed = 77L,
                                  parallel = TRUE, n_cores = 2L)
  bs2 <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, seed = 77L,
                                  parallel = TRUE, n_cores = 2L)

  expect_equal(bs1$adjacency_matrices, bs2$adjacency_matrices)
})

test_that("parallel bootstrap with different seeds gives different results", {
  skip_if_not(parallel::detectCores() >= 2L, "requires >= 2 cores")

  dat <- generate_lingam_sample_6(n = 300, seed = 1)

  bs1 <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, seed = 10L,
                                  parallel = TRUE, n_cores = 2L)
  bs2 <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, seed = 20L,
                                  parallel = TRUE, n_cores = 2L)

  expect_false(identical(bs1$adjacency_matrices, bs2$adjacency_matrices))
})

test_that("parallel and sequential results differ (L'Ecuyer vs set.seed)", {
  # As documented, parallel and sequential runs do not produce identical
  # numbers even with the same seed.
  # This test records that specification as a regression test.
  skip_if_not(parallel::detectCores() >= 2L, "requires >= 2 cores")

  dat <- generate_lingam_sample_6(n = 300, seed = 1)

  bs_seq <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, seed = 42L,
                                    parallel = FALSE)
  bs_par <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, seed = 42L,
                                    parallel = TRUE, n_cores = 2L)

  expect_false(identical(bs_seq$adjacency_matrices, bs_par$adjacency_matrices))
})

test_that("parallel bootstrap returns same structure as sequential", {
  skip_if_not(parallel::detectCores() >= 2L, "requires >= 2 cores")

  dat <- generate_lingam_sample_6(n = 300, seed = 1)

  bs_seq <- lingam_direct_bootstrap(dat$data, n_sampling = 8L, seed = 1L,
                                    parallel = FALSE)
  bs_par <- lingam_direct_bootstrap(dat$data, n_sampling = 8L, seed = 1L,
                                    parallel = TRUE, n_cores = 2L)

  # Even if the numbers differ, the structure (dimensions/class) must be the same
  expect_s3_class(bs_par, "BootstrapResult")
  expect_equal(dim(bs_par$adjacency_matrices), dim(bs_seq$adjacency_matrices))
})
