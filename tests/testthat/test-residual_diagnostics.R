test_that("residual-based diagnostics validate lingam_result", {
  dat <- generate_lingam_sample_6(n = 100, seed = 1)

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
  dat <- generate_lingam_sample_6(n = 100, seed = 1)
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
  dat <- generate_lingam_sample_6(n = 100, seed = 1)
  res <- lingam_direct(dat$data, reg_method = "ols")

  expect_error(get_error_independence_p_values(dat$data, res, method = "bad"))
})
