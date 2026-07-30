# Minimal VARMALiNGAMResult for stationarity checks with known coefficients
# (check_varma_stationarity only reads ar_coefs / ma_coefs / order).
make_varma_result <- function(Phi1, Theta1) {
  structure(
    list(
      ar_coefs = array(Phi1, c(1, 3, 3)),
      ma_coefs = array(Theta1, c(1, 3, 3)),
      order = c(1L, 1L)
    ),
    class = "VARMALiNGAMResult"
  )
}

test_that("check_varma_stationarity reports a stationary invertible fit", {
  st <- check_varma_stationarity(fit_varma_1000())

  expect_s3_class(st, "varma_stationarity")
  expect_true(st$is_stationary)
  expect_true(st$is_invertible)
  expect_equal(st$order, c(1L, 1L))
  expect_length(st$ar_moduli, 3L)
  expect_length(st$ma_moduli, 3L)
  # sorted descending, max fields consistent
  expect_equal(st$ar_moduli, sort(st$ar_moduli, decreasing = TRUE))
  expect_equal(st$max_ar_modulus, st$ar_moduli[1])
  expect_equal(st$max_ma_modulus, st$ma_moduli[1])

  out <- capture.output(print(st))
  expect_true(any(grepl("Stationary:           YES", out)))
  expect_true(any(grepl("Invertible:           YES", out)))
})

test_that("check_varma_stationarity flags non-stationary / non-invertible coefficients", {
  # explosive AR
  st_ar <- check_varma_stationarity(make_varma_result(1.2 * diag(3), 0.2 * diag(3)))
  expect_false(st_ar$is_stationary)
  expect_true(st_ar$is_invertible)
  out <- capture.output(print(st_ar))
  expect_true(any(grepl("non-stationary", out)))

  # non-invertible MA
  st_ma <- check_varma_stationarity(make_varma_result(0.3 * diag(3), 1.1 * diag(3)))
  expect_true(st_ma$is_stationary)
  expect_false(st_ma$is_invertible)
  out <- capture.output(print(st_ma))
  expect_true(any(grepl("not invertible", out)))
})

test_that("check_varma_stationarity handles pure AR models (q = 0)", {
  m <- lingam_varma(varmas_1000_s42()$data,
    order = c(1, 0), criterion = NULL,
    reg_method = "ols", prune = FALSE
  )
  st <- check_varma_stationarity(m)
  expect_length(st$ma_moduli, 0L)
  expect_equal(st$max_ma_modulus, 0)
  expect_true(st$is_invertible)
})

test_that("varma diagnostics validate the result class", {
  expect_error(check_varma_stationarity(list()), "VARMALiNGAMResult")
  expect_error(test_varmalingam_residual_normality(list()), "VARMALiNGAMResult")
  expect_error(plot_varmalingam_residual_qq(list()), "VARMALiNGAMResult")
})

test_that("test_varmalingam_residual_normality rejects normality for uniform errors", {
  m <- fit_varma_1000()
  r <- test_varmalingam_residual_normality(m)

  expect_s3_class(r, "lingam_normality_test")
  expect_equal(nrow(r), 3L)
  # uniform errors: normality should be rejected for every variable
  expect_true(all(r$p_value < 0.05))

  # the innovations e_t and the reduced-form residuals n_t differ
  r_varma <- test_varmalingam_residual_normality(m, on = "varma")
  expect_false(identical(r$p_value, r_varma$p_value))
})

test_that("test_varmalingam_residual_normality_all aggregates methods", {
  m <- fit_varma_1000()
  methods <- "shapiro"
  if (requireNamespace("tseries", quietly = TRUE)) methods <- c("shapiro", "jb")

  out <- test_varmalingam_residual_normality_all(m, methods = methods)
  expect_true(all(c("variable", "skewness", "kurtosis", "p_shapiro", "all_non_gauss")
  %in% names(out)))
  expect_equal(nrow(out), 3L)
  expect_true(all(out$all_non_gauss))
})

test_that("plot_varmalingam_residual_qq returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  m <- fit_varma_1000()
  p <- plot_varmalingam_residual_qq(m)
  expect_s3_class(p, "ggplot")
  p2 <- plot_varmalingam_residual_qq(m, on = "varma")
  expect_s3_class(p2, "ggplot")
})
