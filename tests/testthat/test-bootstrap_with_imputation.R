# NA 入りデータの生成は helper-fixtures.R の make_missing_sample6() に共通化した。

test_that("bootstrap_with_imputation returns ImputationBootstrapResult with expected structure", {
  skip_if_not_installed("mice")
  d <- make_missing_sample6()

  res <- bootstrap_with_imputation(d$X, n_sampling = 5L, n_repeats = 3L, seed = 42, verbose = FALSE)

  expect_s3_class(res, "ImputationBootstrapResult")
  expect_named(res, c(
    "causal_orders", "adjacency_matrices", "resampled_indices",
    "imputation_results", "n_missing"
  ))
  expect_equal(dim(res$causal_orders), c(5L, 6L))
  expect_equal(dim(res$adjacency_matrices), c(5L, 3L, 6L, 6L))
  expect_equal(dim(res$resampled_indices), c(5L, 300L))
  expect_equal(dim(res$imputation_results), c(5L, 3L, 300L, 6L))
  expect_equal(res$n_missing, length(d$na_idx))
})

test_that("imputation_results is non-NA exactly at the bootstrap sample's missing positions", {
  skip_if_not_installed("mice")
  d <- make_missing_sample6()

  res <- bootstrap_with_imputation(d$X, n_sampling = 4L, n_repeats = 2L, seed = 1, verbose = FALSE)

  for (i in seq_len(nrow(res$resampled_indices))) {
    idx <- res$resampled_indices[i, ]
    X_boot <- as.matrix(d$X)[idx, , drop = FALSE]
    expected_non_na_mask <- is.na(X_boot) # non-NA in the result <=> NA in X_boot
    for (r in 1:2) {
      slice <- res$imputation_results[i, r, , ]
      actual_non_na_mask <- !is.na(slice)
      expect_equal(actual_non_na_mask, expected_non_na_mask, ignore_attr = TRUE)
    }
  }
})

test_that("bootstrap_with_imputation is reproducible with the same seed", {
  skip_if_not_installed("mice")
  d <- make_missing_sample6()

  res1 <- bootstrap_with_imputation(d$X, n_sampling = 4L, n_repeats = 2L, seed = 123, verbose = FALSE)
  res2 <- bootstrap_with_imputation(d$X, n_sampling = 4L, n_repeats = 2L, seed = 123, verbose = FALSE)

  expect_equal(res1$causal_orders, res2$causal_orders)
  expect_equal(res1$adjacency_matrices, res2$adjacency_matrices)
  expect_equal(res1$imputation_results, res2$imputation_results)
})

test_that("bootstrap_with_imputation warns when X has no missing values", {
  dat <- sample6_100()$data

  expect_warning(
    bootstrap_with_imputation(dat,
      n_sampling = 2L, n_repeats = 2L,
      imputer = function(X_boot) list(X_boot, X_boot),
      cd_fit = function(X_list) {
        list(causal_order = seq_len(ncol(X_list[[1]])), adjacency_matrices = list(
          matrix(0, ncol(X_list[[1]]), ncol(X_list[[1]])),
          matrix(0, ncol(X_list[[1]]), ncol(X_list[[1]]))
        ))
      },
      verbose = FALSE
    ),
    "no missing values"
  )
})

test_that("custom imputer works, and an invalid one raises a clear error", {
  d <- make_missing_sample6()
  p <- ncol(d$X)

  mean_impute <- function(X_boot) {
    complete <- X_boot
    for (j in seq_len(ncol(complete))) {
      col <- complete[, j]
      col[is.na(col)] <- mean(col, na.rm = TRUE)
      complete[, j] <- col
    }
    list(complete, complete) # n_repeats = 2, ignoring the argument
  }
  fake_cd_fit <- function(X_list) {
    list(
      causal_order = seq_len(p),
      adjacency_matrices = lapply(X_list, function(x) matrix(0, p, p))
    )
  }

  res <- bootstrap_with_imputation(d$X,
    n_sampling = 3L, imputer = mean_impute, cd_fit = fake_cd_fit,
    seed = 1, verbose = FALSE
  )
  expect_s3_class(res, "ImputationBootstrapResult")
  expect_equal(dim(res$adjacency_matrices)[2], 2L)

  # imputer leaves NA behind -> violates its specification
  bad_imputer_na <- function(X_boot) list(X_boot, X_boot) # still has NA
  expect_error(
    bootstrap_with_imputation(d$X,
      n_sampling = 1L, imputer = bad_imputer_na, cd_fit = fake_cd_fit,
      seed = 1, verbose = FALSE
    ),
    "violates its specification"
  )

  # imputer returns wrong dimensions -> violates its specification
  bad_imputer_dim <- function(X_boot) list(X_boot[, 1:2, drop = FALSE])
  expect_error(
    bootstrap_with_imputation(d$X,
      n_sampling = 1L, imputer = bad_imputer_dim, cd_fit = fake_cd_fit,
      seed = 1, verbose = FALSE
    ),
    "violates its specification"
  )
})

test_that("custom cd_fit works, and an invalid one raises a clear error", {
  skip_if_not_installed("mice")
  d <- make_missing_sample6()

  custom_cd_fit <- function(X_list) {
    res <- lingam_multi_group(X_list, reg_method = "ols")
    list(causal_order = res$causal_order, adjacency_matrices = res$adjacency_matrices)
  }

  res <- bootstrap_with_imputation(d$X,
    n_sampling = 2L, n_repeats = 2L, cd_fit = custom_cd_fit,
    seed = 1, verbose = FALSE
  )
  expect_s3_class(res, "ImputationBootstrapResult")

  bad_cd_fit <- function(X_list) {
    list(causal_order = c(1L, 1L, 2L, 3L, 4L, 5L), adjacency_matrices = lapply(X_list, function(x) matrix(0, 6, 6)))
  }
  expect_error(
    bootstrap_with_imputation(d$X,
      n_sampling = 1L, n_repeats = 2L, cd_fit = bad_cd_fit,
      seed = 1, verbose = FALSE
    ),
    "violates its specification"
  )
})

test_that("as_bootstrap_result converts to a usable BootstrapResult", {
  skip_if_not_installed("mice")
  d <- make_missing_sample6()

  res <- bootstrap_with_imputation(d$X, n_sampling = 4L, n_repeats = 3L, seed = 7, verbose = FALSE)
  bs <- as_bootstrap_result(res, aggregate = "median")

  expect_s3_class(bs, "BootstrapResult")
  expect_equal(dim(bs$adjacency_matrices), c(4L, 6L, 6L))
  expect_null(bs$total_effects)
  expect_no_error(get_probabilities(bs))
  expect_error(get_total_causal_effects(bs), "compute_total_effects = FALSE")

  bs_mean <- as_bootstrap_result(res, aggregate = "mean")
  expect_s3_class(bs_mean, "BootstrapResult")
})

test_that("bootstrap_with_imputation validates its arguments", {
  d <- make_missing_sample6()

  expect_error(bootstrap_with_imputation(d$X, n_sampling = 0L), "n_sampling must be a positive integer")
  expect_error(bootstrap_with_imputation(d$X, n_sampling = 2L, n_repeats = 0L), "n_repeats must be a positive integer")
  expect_error(
    bootstrap_with_imputation(d$X, n_sampling = 2L, prior_knowledge = matrix(0, 2, 2)),
    "shape of prior knowledge"
  )
})

test_that("print.ImputationBootstrapResult reports the expected header", {
  skip_if_not_installed("mice")
  d <- make_missing_sample6()

  res <- bootstrap_with_imputation(d$X, n_sampling = 2L, n_repeats = 2L, seed = 1, verbose = FALSE)
  expect_output(print(res), "ImputationBootstrapResult")
})
