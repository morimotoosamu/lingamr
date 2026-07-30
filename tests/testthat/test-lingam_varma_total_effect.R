test_that("estimate_varma_total_effect recovers known effects", {
  s <- varmas_1000_s42()
  m <- fit_varma_1000()

  # contemporaneous: x0 -> x1 -> x2, total = 0.6 * (-0.5) = -0.3
  te <- estimate_varma_total_effect(s$data, m, from_index = 1, to_index = 3)
  expect_lt(abs(te - (-0.3)), 0.1)

  # lag-1 effects: entry [i, j] of Phi1 (I - B0)^{-1} is the total effect of
  # x_j(t-1) on x_i(t) (lag edge plus instantaneous propagation at time t)
  M <- s$true_phi1 %*% solve(diag(3) - s$true_B0)
  te_31 <- estimate_varma_total_effect(s$data, m, from_index = 3, to_index = 1, from_lag = 1)
  expect_lt(abs(te_31 - M[1, 3]), 0.1) # true 0.2
  te_13 <- estimate_varma_total_effect(s$data, m, from_index = 1, to_index = 3, from_lag = 1)
  expect_lt(abs(te_13 - M[3, 1]), 0.1) # true -0.09
})

test_that("estimate_varma_total_effect resolves variable names", {
  s <- varmas_1000_s42()
  m <- fit_varma_1000()

  by_index <- estimate_varma_total_effect(s$data, m, 1, 3)
  by_name <- estimate_varma_total_effect(s$data, m, "x0", "x2")
  expect_identical(by_index, by_name)

  expect_error(
    estimate_varma_total_effect(s$data, m, "nope", "x2"),
    "not found"
  )
})

test_that("estimate_varma_total_effect warns on reversed causal order", {
  s <- varmas_1000_s42()
  m <- fit_varma_1000()

  # x2 is last in the causal order, so x2 -> x0 at lag 0 is reversed
  expect_warning(
    estimate_varma_total_effect(s$data, m, from_index = 3, to_index = 1),
    "Causal order"
  )
  # no warning when the source is lagged
  expect_silent(
    estimate_varma_total_effect(s$data, m, from_index = 3, to_index = 1, from_lag = 1)
  )
})

test_that("estimate_varma_total_effect validates inputs", {
  s <- varmas_1000_s42()
  m <- fit_varma_1000()

  expect_error(
    estimate_varma_total_effect(s$data, list(), 1, 3),
    "VARMALiNGAMResult"
  )
  expect_error(
    estimate_varma_total_effect(s$data[, 1:2], m, 1, 2),
    "2 variables but result"
  )
  expect_error(
    estimate_varma_total_effect(s$data[1:100, ], m, 1, 3),
    "row count"
  )
  expect_error(
    estimate_varma_total_effect(s$data, m, 1, 3, from_lag = -1),
    "from_lag"
  )
  expect_error(
    estimate_varma_total_effect(s$data, m, 5, 3),
    "between 1 and 3"
  )
})
