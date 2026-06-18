test_that("lingam_var returns VARLiNGAMResult with correct structure", {
  s <- generate_varlingam_sample(n = 500, seed = 42)
  m <- lingam_var(s$data, lags = 1, reg_method = "ols", criterion = NULL, prune = FALSE)

  expect_s3_class(m, "VARLiNGAMResult")
  expect_named(m, c("adjacency_matrices", "causal_order", "residuals", "lags"))
  expect_equal(dim(m$adjacency_matrices), c(2L, 3L, 3L))  # 1 + lags slices
  expect_equal(m$lags, 1L)
  expect_equal(length(m$causal_order), 3L)
})

test_that("lingam_var recovers the instantaneous structure (x0 -> x1 -> x2)", {
  s <- generate_varlingam_sample(n = 2000, seed = 42)
  m <- lingam_var(s$data, lags = 1, reg_method = "ols", criterion = NULL, prune = FALSE)
  B0 <- m$adjacency_matrices[1, , ]

  # x1 <- x0 is positive (true 0.6), x2 <- x1 is negative (true -0.5)
  expect_gt(B0[2, 1], 0.3)
  expect_lt(B0[3, 2], -0.3)
  expect_equal(m$causal_order, c(1L, 2L, 3L))
})

test_that("lingam_var lag-1 matrix recovers the true M1", {
  s <- generate_varlingam_sample(n = 2000, seed = 42)
  m <- lingam_var(s$data, lags = 1, reg_method = "ols", criterion = NULL, prune = FALSE)
  B1 <- m$adjacency_matrices[2, , ]

  # diagonal of B1 should be close to the true M1 diagonal (0.4, 0.3, 0.5)
  expect_equal(unname(diag(B1)), c(0.4, 0.3, 0.5), tolerance = 0.1)
})

test_that("lingam_var selects the lag order by BIC", {
  s <- generate_varlingam_sample(n = 2000, seed = 42)
  m <- lingam_var(s$data, lags = 5, reg_method = "ols", criterion = "bic", prune = FALSE)

  expect_equal(m$lags, 1L)  # true lag order is 1
})

test_that("lingam_var supports aic / hqic / fpe lag selection", {
  s <- generate_varlingam_sample(n = 2000, seed = 42)
  # All criteria should recover the true lag order (1) on this lag-1 model.
  for (crit in c("aic", "hqic", "fpe")) {
    m <- lingam_var(s$data, lags = 5, reg_method = "ols", criterion = crit, prune = FALSE)
    expect_equal(m$lags, 1L, info = crit)
  }
})

test_that("select_var_lag compares candidates on a common sample", {
  # On a common sample the residual covariance of each candidate is computed
  # over t = max_lag + 1 .. n, so the selected lag is invariant to max_lag
  # as long as the true lag (1) stays within range.
  s <- generate_varlingam_sample(n = 2000, seed = 42)
  X <- as.matrix(s$data)
  expect_equal(select_var_lag(X, max_lag = 3, criterion = "bic"), 1L)
  expect_equal(select_var_lag(X, max_lag = 8, criterion = "bic"), 1L)
})

test_that("lingam_var errors on missing values", {
  s <- generate_varlingam_sample(n = 200, seed = 1)
  X <- s$data
  X[5, 2] <- NA
  expect_error(lingam_var(X, lags = 1, reg_method = "ols", criterion = NULL, prune = FALSE), "NA")
})

test_that("lingam_var errors on invalid inputs", {
  s <- generate_varlingam_sample(n = 200, seed = 1)

  expect_error(lingam_var(s$data[, 1, drop = FALSE], reg_method = "ols"))   # 1 variable
  expect_error(lingam_var(s$data, lags = 0, reg_method = "ols"))            # lags < 1
  expect_error(lingam_var(s$data, lags = 1, criterion = "bad", reg_method = "ols"))
})

test_that("print.VARLiNGAMResult runs without error", {
  s <- generate_varlingam_sample(n = 300, seed = 1)
  m <- lingam_var(s$data, lags = 1, reg_method = "ols", criterion = NULL, prune = FALSE)

  expect_output(print(m), "VAR-LiNGAM Result")
  expect_output(print(m), "Lag order")
  expect_output(print(m), "Causal order")
})

test_that("generate_varlingam_sample returns the expected structure", {
  s <- generate_varlingam_sample(n = 100, seed = 1)

  expect_named(s, c("data", "true_B0", "true_M1"))
  expect_equal(nrow(s$data), 100L)
  expect_equal(ncol(s$data), 3L)
  expect_equal(dim(s$true_B0), c(3L, 3L))
  expect_equal(dim(s$true_M1), c(3L, 3L))
})

test_that("lingam_var validates the prune argument", {
  s <- generate_varlingam_sample(n = 100, seed = 1)
  expect_error(
    lingam_var(s$data, lags = 1, reg_method = "ols", criterion = NULL, prune = "yes"),
    "prune"
  )
})

test_that("prune = TRUE recovers the structure and shrinks weak edges", {
  skip_if_not_installed("glmnet")
  s <- generate_varlingam_sample(n = 2000, seed = 42)
  m <- lingam_var(s$data, lags = 1, reg_method = "ols", criterion = NULL, prune = TRUE)
  B0 <- m$adjacency_matrices[1, , ]
  B1 <- m$adjacency_matrices[2, , ]

  # True instantaneous edges survive pruning.
  expect_gt(B0[2, 1], 0.3)   # x0 -> x1 (true 0.6)
  expect_lt(B0[3, 2], -0.3)  # x1 -> x2 (true -0.5)
  # Lag-1 diagonal (true 0.4, 0.3, 0.5) is retained.
  expect_equal(unname(diag(B1)), c(0.4, 0.3, 0.5), tolerance = 0.15)
  # A structurally-absent instantaneous edge (x2 -> x0) is pruned to zero.
  expect_equal(B0[1, 3], 0)
})
