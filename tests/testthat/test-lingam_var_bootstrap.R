# Fast settings throughout: OLS instantaneous structure, no pruning, no glmnet.

test_that("lingam_var_bootstrap returns a well-formed VARBootstrapResult", {
  s <- generate_varlingam_sample(n = 500, seed = 42)
  bs <- lingam_var_bootstrap(s$data,
    n_sampling = 20L, lags = 1, criterion = NULL,
    reg_method = "ols", prune = FALSE, seed = 1, verbose = FALSE
  )

  expect_s3_class(bs, "VARBootstrapResult")
  expect_equal(bs$lags, 1L)
  # adjacency: list of length n_sampling, each n_features x n_features*(1 + lags)
  expect_length(bs$adjacency_matrices, 20L)
  expect_equal(dim(bs$adjacency_matrices[[1]]), c(3L, 6L))
  # total effects array
  expect_equal(dim(bs$total_effects), c(20L, 3L, 6L))
  expect_equal(dim(bs$causal_orders), c(20L, 3L))
  expect_output(print(bs), "VARBootstrapResult")
})

test_that("lingam_var_bootstrap is reproducible for a fixed seed (sequential)", {
  s <- generate_varlingam_sample(n = 400, seed = 7)
  bs1 <- lingam_var_bootstrap(s$data,
    n_sampling = 10L, criterion = NULL,
    reg_method = "ols", prune = FALSE, seed = 123, verbose = FALSE
  )
  bs2 <- lingam_var_bootstrap(s$data,
    n_sampling = 10L, criterion = NULL,
    reg_method = "ols", prune = FALSE, seed = 123, verbose = FALSE
  )
  expect_equal(bs1$adjacency_matrices[[1]], bs2$adjacency_matrices[[1]])
  expect_equal(bs1$total_effects, bs2$total_effects)
})

test_that("get_var_probabilities has the right shape and detects strong edges", {
  s <- generate_varlingam_sample(n = 800, seed = 42)
  bs <- lingam_var_bootstrap(s$data,
    n_sampling = 20L, criterion = NULL,
    reg_method = "ols", prune = FALSE, seed = 1, verbose = FALSE
  )

  p <- get_var_probabilities(bs, min_causal_effect = 0.3)
  expect_equal(dim(p), c(3L, 6L))
  expect_true(all(p >= 0 & p <= 1))
  # x0 -> x1 (true 0.6) is a strong, frequently-detected contemporaneous edge
  expect_gt(p[2, 1], 0.8)
  # x2 -> x0 contemporaneously is structurally absent (root has no parents)
  p0 <- get_var_probabilities(bs)  # threshold 0
  expect_equal(p0[1, 3], 0)
})

test_that("get_var_paths finds the indirect path x0 -> x1 -> x2", {
  s <- generate_varlingam_sample(n = 800, seed = 42)
  bs <- lingam_var_bootstrap(s$data,
    n_sampling = 20L, criterion = NULL,
    reg_method = "ols", prune = FALSE, seed = 1, verbose = FALSE
  )

  paths <- get_var_paths(bs, from_index = 1, to_index = 3)
  expect_s3_class(paths, "data.frame")
  expect_named(paths, c("path", "effect", "probability"))
  expect_gt(nrow(paths), 0)
  expect_true(all(paths$probability > 0 & paths$probability <= 1))
  # the contemporaneous chain 1 -> 2 -> 3 should be among the enumerated paths
  has_chain <- any(vapply(paths$path, function(p) identical(p, c(1L, 2L, 3L)), logical(1)))
  expect_true(has_chain)
})

test_that("bootstrap and path helpers validate inputs", {
  s <- generate_varlingam_sample(n = 300, seed = 1)
  expect_error(
    lingam_var_bootstrap(s$data, n_sampling = 0L, criterion = NULL,
      reg_method = "ols", prune = FALSE, verbose = FALSE),
    "n_sampling"
  )
  bs <- lingam_var_bootstrap(s$data,
    n_sampling = 5L, criterion = NULL,
    reg_method = "ols", prune = FALSE, seed = 1, verbose = FALSE
  )
  # to_lag must not exceed from_lag
  expect_error(get_var_paths(bs, 1, 2, from_lag = 0, to_lag = 1), "from_lag")
  # same variable, same lag is invalid
  expect_error(get_var_paths(bs, 1, 1, from_lag = 0, to_lag = 0), "same")
})

test_that("lingam_var_bootstrap runs with pruning (glmnet)", {
  skip_if_not_installed("glmnet")
  s <- generate_varlingam_sample(n = 400, seed = 42)
  bs <- lingam_var_bootstrap(s$data,
    n_sampling = 5L, criterion = NULL,
    reg_method = "adaptive_lasso", prune = TRUE, seed = 1, verbose = FALSE
  )
  expect_s3_class(bs, "VARBootstrapResult")
  expect_equal(dim(bs$total_effects), c(5L, 3L, 6L))
})
