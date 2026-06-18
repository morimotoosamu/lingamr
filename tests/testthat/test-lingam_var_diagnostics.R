# Helper: build a minimal VARLiNGAMResult (lag 1) for targeted stationarity tests.
make_var_result <- function(B0, B1) {
  p <- ncol(B0)
  am <- array(0, dim = c(2L, p, p))
  am[1, , ] <- B0
  am[2, , ] <- B1
  structure(
    list(
      adjacency_matrices = am,
      causal_order = seq_len(p),
      residuals = matrix(0, nrow = 10L, ncol = p),
      lags = 1L
    ),
    class = "VARLiNGAMResult"
  )
}

test_that("check_var_stationarity flags a stationary fitted model", {
  s <- generate_varlingam_sample(n = 1000, seed = 42)
  m <- lingam_var(s$data, lags = 1, reg_method = "ols", criterion = NULL, prune = FALSE)
  st <- check_var_stationarity(m)

  expect_s3_class(st, "var_stationarity")
  expect_true(st$is_stationary)
  expect_lt(st$max_modulus, 1)
  expect_length(st$moduli, 3L)  # p * lags = 3
  expect_output(print(st), "Stationarity Check")
})

test_that("check_var_stationarity detects a non-stationary model", {
  p <- 3L
  B0 <- matrix(0, p, p)
  # explosive lag-1 dynamics: a diagonal entry > 1 puts a root outside the circle
  B1 <- diag(c(1.2, 0.5, 0.5))
  res <- make_var_result(B0, B1)
  st <- check_var_stationarity(res)

  expect_false(st$is_stationary)
  expect_gt(st$max_modulus, 1)
  expect_output(print(st), "non-stationary")
})

test_that("check_var_stationarity validates its input", {
  expect_error(check_var_stationarity(list()), "VARLiNGAMResult")
})

test_that("test_varlingam_residual_normality rejects normality for uniform errors", {
  # generate_varlingam_sample uses uniform (non-Gaussian) errors
  s <- generate_varlingam_sample(n = 1500, seed = 42)
  m <- lingam_var(s$data, lags = 1, reg_method = "ols", criterion = NULL, prune = FALSE)

  nt <- test_varlingam_residual_normality(m, method = "shapiro", on = "innovations")
  expect_s3_class(nt, "lingam_normality_test")
  expect_equal(nrow(nt), 3L)
  # uniform innovations are strongly non-Gaussian (negative excess kurtosis)
  expect_true(all(nt$is_non_gauss))
  expect_true(all(nt$kurtosis < 0))
  expect_output(print(nt), "Residual Normality Test")
})

test_that("test_varlingam_residual_normality supports both targets", {
  s <- generate_varlingam_sample(n = 1000, seed = 1)
  m <- lingam_var(s$data, lags = 1, reg_method = "ols", criterion = NULL, prune = FALSE)

  innov <- test_varlingam_residual_normality(m, on = "innovations")
  varr <- test_varlingam_residual_normality(m, on = "var")
  expect_s3_class(innov, "lingam_normality_test")
  expect_s3_class(varr, "lingam_normality_test")
  # the two targets generally differ (innovations are (I - B0) n_t)
  expect_false(isTRUE(all.equal(innov$kurtosis, varr$kurtosis)))
})

test_that("test_varlingam_residual_normality validates its input", {
  expect_error(test_varlingam_residual_normality(list()), "VARLiNGAMResult")
})

test_that("test_varlingam_residual_normality_all builds a multi-method table", {
  skip_if_not_installed("nortest")
  skip_if_not_installed("tseries")
  s <- generate_varlingam_sample(n = 1500, seed = 42)
  m <- lingam_var(s$data, lags = 1, reg_method = "ols", criterion = NULL, prune = FALSE)

  tbl <- test_varlingam_residual_normality_all(m, methods = c("shapiro", "ad", "lillie", "jb"))
  expect_s3_class(tbl, "data.frame")
  expect_equal(nrow(tbl), 3L)
  # one p-value column per requested method, plus base stats and summary
  expect_true(all(c("p_shapiro", "p_ad", "p_lillie", "p_jb") %in% names(tbl)))
  expect_true(all(c("variable", "skewness", "kurtosis", "all_non_gauss") %in% names(tbl)))
  # uniform innovations are non-Gaussian under every test
  expect_true(all(tbl$all_non_gauss))
})

test_that("test_varlingam_residual_normality_all validates input", {
  expect_error(test_varlingam_residual_normality_all(list()), "VARLiNGAMResult")
})

test_that("plot_varlingam_residual_qq returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  s <- generate_varlingam_sample(n = 500, seed = 1)
  m <- lingam_var(s$data, lags = 1, reg_method = "ols", criterion = NULL, prune = FALSE)

  p_innov <- plot_varlingam_residual_qq(m, on = "innovations")
  p_var <- plot_varlingam_residual_qq(m, on = "var")
  expect_s3_class(p_innov, "ggplot")
  expect_s3_class(p_var, "ggplot")
})

test_that("plot_varlingam_residual_qq validates input", {
  expect_error(plot_varlingam_residual_qq(list()), "VARLiNGAMResult")
})
