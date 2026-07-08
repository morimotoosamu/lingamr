test_that("estimate_total_effect errors on invalid lingam_result", {
  dat <- generate_lingam_sample_6(n = 200, seed = 1)
  fake <- list(adjacency_matrix = matrix(0, 6, 6), causal_order = 1:6)

  expect_error(estimate_total_effect(dat$data, fake, 1, 2), "must be the return value of lingam_direct")
})

test_that("estimate_total_effect errors when X dimensions mismatch", {
  dat <- generate_lingam_sample_6(n = 200, seed = 1)
  res <- lingam_direct(dat$data, reg_method = "ols")
  X_bad <- dat$data[, 1:4]   # different number of variables

  expect_error(estimate_total_effect(X_bad, res, 1, 2), "variables but lingam_result")
})

test_that("estimate_total_effect errors when from_index == to_index", {
  dat <- generate_lingam_sample_6(n = 200, seed = 1)
  res <- lingam_direct(dat$data, reg_method = "ols")

  expect_error(estimate_total_effect(dat$data, res, 1, 1), "must differ")
})

test_that("estimate_total_effect recovers the known direct and indirect effects", {
  # Known DGP (see generate_lingam_sample_6): x3 -> x0 (3.0) -> x1 (3.0),
  #                                            x3 -> x2 (6.0) -> x1 (2.0)
  # so the total effect x3 -> x1 is the sum over both paths: 3*3 + 6*2 = 21.
  dat <- generate_lingam_sample_6(n = 2000, seed = 42)
  res <- lingam_direct(dat$data, reg_method = "ols")

  te_direct <- estimate_total_effect(dat$data, res, "x3", "x0", method = "ols")
  expect_equal(unname(te_direct), 3.0, tolerance = 0.2)

  te_indirect <- estimate_total_effect(dat$data, res, "x3", "x1", method = "ols")
  expect_equal(unname(te_indirect), 21.0, tolerance = 1.0)
})

test_that("estimate_total_effect accepts variable names", {
  dat <- generate_lingam_sample_6(n = 500, seed = 42)
  res <- lingam_direct(dat$data, reg_method = "ols")

  te_idx  <- estimate_total_effect(dat$data, res, 4, 1)   # x3(4) -> x0(1)
  te_name <- estimate_total_effect(dat$data, res, "x3", "x0")

  expect_equal(te_idx, te_name)
})

test_that("estimate_all_total_effects returns correctly shaped matrix", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  res <- lingam_direct(dat$data, reg_method = "ols")
  TE  <- estimate_all_total_effects(dat$data, res, method = "ols")

  expect_true(is.matrix(TE))
  expect_equal(dim(TE), c(6L, 6L))
  # diagonal is zero (no self-effect)
  expect_true(all(diag(TE) == 0))
  # x3 is exogenous, so its column (cause side) should be zero
  expect_true(all(TE["x3", ] == 0))
})

test_that("estimate_all_total_effects recovers known effect sizes", {
  dat <- generate_lingam_sample_6(n = 2000, seed = 42)
  res <- lingam_direct(dat$data, reg_method = "ols")
  TE  <- estimate_all_total_effects(dat$data, res, method = "ols")

  expect_equal(TE["x0", "x3"], 3.0, tolerance = 0.2)   # x3 -> x0, direct
  expect_equal(TE["x1", "x3"], 21.0, tolerance = 1.0)  # x3 -> x1, via x0 and x2
  expect_equal(TE["x4", "x0"], 8.0, tolerance = 0.5)   # x0 -> x4, direct
})

test_that("estimate_all_total_effects errors on invalid lingam_result", {
  dat <- generate_lingam_sample_6(n = 200, seed = 1)
  fake <- list(adjacency_matrix = matrix(0, 6, 6), causal_order = 1:6)

  expect_error(estimate_all_total_effects(dat$data, fake), "must be the return value of lingam_direct")
})
