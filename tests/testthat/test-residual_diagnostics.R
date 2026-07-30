test_that("residual-based diagnostics validate lingam_result", {
  dat <- sample6_100()

  expect_error(
    get_error_independence_p_values(dat$data, "not a result"),
    "lingam_direct"
  )
  expect_error(
    test_residual_normality(dat$data, list(adjacency_matrix = diag(6))),
    "lingam_direct"
  )
})

test_that("residual-based diagnostics validate dimension mismatch", {
  dat <- sample6_100()
  res <- lingam_direct(dat$data, reg_method = "ols")

  expect_error(
    get_error_independence_p_values(dat$data[, 1:5], res),
    "variables"
  )
  expect_error(
    test_residual_normality(dat$data[, 1:5], res),
    "variables"
  )
})

test_that("get_error_independence_p_values validates method", {
  dat <- sample6_100()
  res <- lingam_direct(dat$data, reg_method = "ols")

  expect_error(get_error_independence_p_values(dat$data, res, method = "bad"))
})

test_that("get_error_independence_p_values returns a valid symmetric p-value matrix", {
  dat <- sample6_300()
  res <- fit_direct_300()

  p <- get_error_independence_p_values(dat$data, res)
  expect_true(is.matrix(p))
  expect_equal(dim(p), c(6L, 6L))
  expect_true(all(diag(p) |> is.na()))
  off_diag <- p[upper.tri(p)]
  expect_true(all(off_diag >= 0 & off_diag <= 1))
  expect_equal(p, t(p))
})

test_that("test_residual_normality detects non-Gaussianity for uniform noise", {
  dat <- generate_lingam_sample_6(n = 1000, seed = 1, noise_dist = "uniform")
  res <- lingam_direct(dat$data, reg_method = "ols")

  out <- test_residual_normality(dat$data, res)
  expect_s3_class(out, "lingam_normality_test")
  expect_true(any(out$is_non_gauss))
})

test_that("plot_residual_qq returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  dat <- sample6_300()
  res <- fit_direct_300()

  p <- plot_residual_qq(dat$data, res)
  expect_s3_class(p, "ggplot")
})

test_that("shapiro subsampling for n > 5000 is deterministic and leaves RNG untouched", {
  dat <- generate_lingam_sample_6(n = 6000, seed = 1)
  res <- lingam_direct(dat$data, reg_method = "ols")

  set.seed(123)
  rs_before <- .Random.seed
  out1 <- suppressWarnings(test_residual_normality(dat$data, res))
  rs_after <- .Random.seed
  out2 <- suppressWarnings(test_residual_normality(dat$data, res))

  # identical p-values across calls, and the caller's RNG stream is untouched
  expect_identical(out1$p_value, out2$p_value)
  expect_identical(rs_before, rs_after)
  # the warning explains the deterministic subsampling
  expect_warning(test_residual_normality(dat$data, res), "deterministic")
})
