# File-local fixture: one light bootstrap shared by the tests below
# (deterministic: seeded data + seeded resampling).
varma_bs_data <- generate_varmalingam_sample(n = 500, seed = 42)
varma_bs <- lingam_varma_bootstrap(varma_bs_data$data,
  n_sampling = 10L, order = c(1, 1), criterion = NULL,
  reg_method = "ols", prune = FALSE, seed = 1, verbose = FALSE
)

test_that("lingam_varma_bootstrap returns the documented structure", {
  bs <- varma_bs
  expect_s3_class(bs, "VARMABootstrapResult")
  expect_named(bs, c(
    "adjacency_matrices", "total_effects", "order",
    "resampled_indices", "causal_orders"
  ))
  expect_length(bs$adjacency_matrices, 10L)
  # joined matrix: 3 x 3*(1 + p + q) = 3 x 9
  expect_equal(dim(bs$adjacency_matrices[[1]]), c(3L, 9L))
  # total effects cover the instantaneous and AR blocks only: 3 x 3*(1 + p)
  expect_equal(dim(bs$total_effects), c(10L, 3L, 6L))
  expect_equal(bs$order, c(1L, 1L))
  expect_equal(dim(bs$causal_orders), c(10L, 3L))
  expect_length(bs$resampled_indices, 10L)
  expect_true(all(vapply(bs$resampled_indices, length, integer(1)) == 500L))

  out <- capture.output(print(bs))
  expect_true(any(grepl("VARMABootstrapResult: 10 samplings, 3 features, order \\(1, 1\\)", out)))
})

test_that("lingam_varma_bootstrap is reproducible with a seed", {
  bs2 <- lingam_varma_bootstrap(varma_bs_data$data,
    n_sampling = 10L, order = c(1, 1), criterion = NULL,
    reg_method = "ols", prune = FALSE, seed = 1, verbose = FALSE
  )
  expect_identical(varma_bs$adjacency_matrices, bs2$adjacency_matrices)
  expect_identical(varma_bs$total_effects, bs2$total_effects)
})

test_that("lingam_varma_bootstrap validates inputs", {
  X <- varma_bs_data$data
  expect_error(lingam_varma_bootstrap(X, n_sampling = 0), "n_sampling")
  expect_error(lingam_varma_bootstrap(X, n_sampling = 5, order = c(0, 0)), "order must be")
  expect_error(
    lingam_varma_bootstrap(X, n_sampling = 5, reg_method = "nope"),
    "should be one of"
  )
  Xna <- X
  Xna[1, 1] <- NA
  expect_error(lingam_varma_bootstrap(Xna, n_sampling = 5), "missing values")
})

test_that("bootstrapped instantaneous effects recover the true structure", {
  # true contemporaneous total effect x0 -> x2 = -0.3
  te_med <- apply(varma_bs$total_effects, c(2, 3), stats::median)
  expect_lt(abs(te_med[3, 1] - (-0.3)), 0.15)
  # true direct edge x0 -> x1 = 0.6
  expect_lt(abs(te_med[2, 1] - 0.6), 0.15)
})

test_that("get_varma_probabilities returns valid probabilities", {
  probs <- get_varma_probabilities(varma_bs)
  expect_equal(dim(probs), c(3L, 9L))
  expect_true(all(probs >= 0 & probs <= 1))

  # strong true edges appear in (nearly) every unpruned OLS sample
  expect_gt(probs[2, 1], 0.9) # x0 -> x1 (instantaneous)
  expect_gt(probs[3, 2], 0.9) # x1 -> x2 (instantaneous)

  # thresholding reduces the probabilities monotonically
  probs_thr <- get_varma_probabilities(varma_bs, min_causal_effect = 0.4)
  expect_true(all(probs_thr <= probs))

  expect_error(get_varma_probabilities(varma_bs, min_causal_effect = -1), ">= 0")
  expect_error(get_varma_probabilities(list()), "VARMABootstrapResult")
})

test_that("get_varma_paths enumerates paths over the psi blocks", {
  paths <- get_varma_paths(varma_bs, from_index = 1, to_index = 3)
  expect_s3_class(paths, "data.frame")
  expect_named(paths, c("path", "effect", "probability"))
  expect_true(all(paths$probability > 0 & paths$probability <= 1))
  # the direct instantaneous chain x0 -> x1 -> x2 should be found in every sample
  path_strs <- vapply(paths$path, paste, "", collapse = "_")
  expect_true("1_2_3" %in% path_strs)
  expect_equal(paths$probability[path_strs == "1_2_3"], 1)

  # lagged source works and node indices shift into the lag block
  paths_lag <- get_varma_paths(varma_bs, from_index = 1, to_index = 3, from_lag = 1)
  expect_true(nrow(paths_lag) > 0)
  expect_true(all(vapply(paths_lag$path, function(p) p[1] == 4L, logical(1))))
})

test_that("get_varma_paths validates lags", {
  expect_error(
    get_varma_paths(varma_bs, 1, 3, from_lag = 2),
    "must not exceed the model's AR order"
  )
  expect_error(get_varma_paths(varma_bs, 1, 3, from_lag = 0, to_lag = 1), "from_lag must be >=")
  expect_error(get_varma_paths(varma_bs, 1, 1), "same variable")
  expect_error(get_varma_paths(varma_bs, 1, 3, min_causal_effect = -1), ">= 0")
})

test_that("lingam_varma_bootstrap works with pruning (glmnet)", {
  skip_if_not_installed("glmnet")
  bs <- lingam_varma_bootstrap(varma_bs_data$data,
    n_sampling = 3L, order = c(1, 1), criterion = NULL,
    reg_method = "ols", prune = TRUE, seed = 1, verbose = FALSE
  )
  expect_length(bs$adjacency_matrices, 3L)
  # pruning zeroes out weak edges, so some entries must be exactly zero
  expect_gt(sum(bs$adjacency_matrices[[1]] == 0), 0L)
})
