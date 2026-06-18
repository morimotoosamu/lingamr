# True model (generate_varlingam_sample):
#   instantaneous: x0 -> x1 (0.6), x1 -> x2 (-0.5)
#   so the total contemporaneous effect x0 -> x2 is 0.6 * -0.5 = -0.3

test_that("estimate_var_total_effect recovers direct and indirect effects", {
  s <- generate_varlingam_sample(n = 3000, seed = 42)
  m <- lingam_var(s$data, lags = 1, reg_method = "ols", criterion = NULL, prune = FALSE)

  # direct contemporaneous effects
  expect_equal(estimate_var_total_effect(s$data, m, 1, 2), 0.6, tolerance = 0.1)
  expect_equal(estimate_var_total_effect(s$data, m, 2, 3), -0.5, tolerance = 0.1)
  # indirect: x0 -> x1 -> x2
  expect_equal(estimate_var_total_effect(s$data, m, 1, 3), -0.3, tolerance = 0.1)
})

test_that("estimate_var_total_effect accepts variable names", {
  s <- generate_varlingam_sample(n = 2000, seed = 42)
  m <- lingam_var(s$data, lags = 1, reg_method = "ols", criterion = NULL, prune = FALSE)

  by_idx <- estimate_var_total_effect(s$data, m, 1, 2)
  by_name <- estimate_var_total_effect(s$data, m, "x0", "x1")
  expect_equal(by_idx, by_name)
})

test_that("estimate_var_total_effect warns on reversed causal order", {
  s <- generate_varlingam_sample(n = 1000, seed = 42)
  m <- lingam_var(s$data, lags = 1, reg_method = "ols", criterion = NULL, prune = FALSE)

  # from = x2, to = x0 contemporaneously is against the causal order
  expect_warning(estimate_var_total_effect(s$data, m, 3, 1), "[Cc]ausal order")
})

test_that("estimate_var_total_effect validates inputs", {
  s <- generate_varlingam_sample(n = 500, seed = 1)
  m <- lingam_var(s$data, lags = 1, reg_method = "ols", criterion = NULL, prune = FALSE)

  expect_error(estimate_var_total_effect(s$data, list(), 1, 2), "VARLiNGAMResult")
  expect_error(estimate_var_total_effect(s$data, m, 1, 2, from_lag = -1), "from_lag")
  expect_error(estimate_var_total_effect(s$data, m, 99, 2), "from_index")
})

test_that("var_total_effect_graph matches the path product", {
  s <- generate_varlingam_sample(n = 2000, seed = 42)
  m <- lingam_var(s$data, lags = 1, reg_method = "ols", criterion = NULL, prune = FALSE)

  # joined adjacency cbind(B0, B1); graph-based x0 -> x2 = product over x0->x1->x2
  am <- m$adjacency_matrices
  am_joined <- cbind(am[1, , ], am[2, , ])
  graph_te <- var_total_effect_graph(am_joined, 1, 3)
  expect_equal(graph_te, -0.3, tolerance = 0.1)
})
