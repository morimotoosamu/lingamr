test_that("estimate_total_effect errors on invalid lingam_result", {
  dat <- generate_lingam_sample_6(n = 200, seed = 1)
  fake <- list(adjacency_matrix = matrix(0, 6, 6), causal_order = 1:6)

  expect_error(estimate_total_effect(dat$data, fake, 1, 2))
})

test_that("estimate_total_effect errors when X dimensions mismatch", {
  dat <- generate_lingam_sample_6(n = 200, seed = 1)
  res <- lingam_direct(dat$data, reg_method = "ols")
  X_bad <- dat$data[, 1:4]   # 変数数が違う

  expect_error(estimate_total_effect(X_bad, res, 1, 2))
})

test_that("estimate_total_effect errors when from_index == to_index", {
  dat <- generate_lingam_sample_6(n = 200, seed = 1)
  res <- lingam_direct(dat$data, reg_method = "ols")

  expect_error(estimate_total_effect(dat$data, res, 1, 1))
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
  # 対角はゼロ（自己効果なし）
  expect_true(all(diag(TE) == 0))
  # x3 は外生変数なので列（原因側）はゼロのはず
  expect_true(all(TE["x3", ] == 0))
})

test_that("estimate_all_total_effects errors on invalid lingam_result", {
  dat <- generate_lingam_sample_6(n = 200, seed = 1)
  fake <- list(adjacency_matrix = matrix(0, 6, 6), causal_order = 1:6)

  expect_error(estimate_all_total_effects(dat$data, fake))
})
