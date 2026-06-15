test_that("summary_lingam returns lingam_summary with expected structure", {
  dat <- generate_lingam_sample_6(n = 500, seed = 42)
  res <- lingam_direct(dat$data, reg_method = "ols")
  s   <- summary_lingam(dat$data, res)

  expect_s3_class(s, "lingam_summary")
  expect_equal(s$n_variables, 6L)
  expect_equal(s$n_samples, 500L)
  expect_equal(length(s$causal_order), 6L)
  expect_true(is.matrix(s$independence_p_values))
  expect_equal(dim(s$independence_p_values), c(6L, 6L))
  # 6 variables -> 15 upper-triangular pairs
  expect_equal(s$n_pairs, 15L)
  expect_true(s$n_dependent_pairs >= 0 && s$n_dependent_pairs <= 15L)
  expect_true(s$n_non_gaussian >= 0 && s$n_non_gaussian <= 6L)
  expect_s3_class(s$normality, "lingam_normality_test")
})

test_that("summary_lingam detects non-Gaussianity for uniform noise", {
  dat <- generate_lingam_sample_6(n = 1000, seed = 1, noise_dist = "uniform")
  res <- lingam_direct(dat$data, reg_method = "ols")
  s   <- summary_lingam(dat$data, res)

  # residuals from a uniform distribution should be detected as non-Gaussian
  expect_gt(s$n_non_gaussian, 0)
})

test_that("print.lingam_summary runs and shows both assumptions", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  res <- lingam_direct(dat$data, reg_method = "ols")
  s   <- summary_lingam(dat$data, res)

  expect_output(print(s), "Direct LiNGAM Model Summary")
  expect_output(print(s), "Assumption 1")
  expect_output(print(s), "Assumption 2")
})

test_that("summary_lingam errors on invalid lingam_result", {
  dat  <- generate_lingam_sample_6(n = 200, seed = 1)
  fake <- list(adjacency_matrix = matrix(0, 6, 6), causal_order = 1:6)

  expect_error(summary_lingam(dat$data, fake))
})

test_that("summary_lingam errors when X dimensions mismatch", {
  dat <- generate_lingam_sample_6(n = 200, seed = 1)
  res <- lingam_direct(dat$data, reg_method = "ols")

  expect_error(summary_lingam(dat$data[, 1:4], res))
})

test_that("summary_lingam validates alpha and method arguments", {
  dat <- generate_lingam_sample_6(n = 200, seed = 1)
  res <- lingam_direct(dat$data, reg_method = "ols")

  expect_error(summary_lingam(dat$data, res, alpha = 0))
  expect_error(summary_lingam(dat$data, res, alpha = 1.5))
  expect_error(summary_lingam(dat$data, res, independence_method = "bad"))
  expect_error(summary_lingam(dat$data, res, normality_method = "bad"))
})
