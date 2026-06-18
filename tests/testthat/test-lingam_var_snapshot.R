# Snapshot (golden-value) regression test for lingam_var().
#
# These pinned numbers are NOT a Python-parity check; they are this R version's
# own deterministic output, frozen to catch *unintended* numerical changes in
# the VAR-LiNGAM pipeline. The configuration below is fully deterministic
# (uniform-seeded data, pwling measure, OLS regression, no pruning, no lag
# selection), so the result is reproducible across runs.
#
# If the algorithm is changed on purpose, regenerate the expected values:
#   s <- generate_varlingam_sample(n = 500, seed = 42)
#   m <- lingam_var(s$data, lags = 1, reg_method = "ols",
#                   criterion = NULL, prune = FALSE)
#   dput(round(unname(m$adjacency_matrices[1, , ]), 7))  # B0
#   dput(round(unname(m$adjacency_matrices[2, , ]), 7))  # B1

test_that("lingam_var output matches the pinned snapshot", {
  s <- generate_varlingam_sample(n = 500, seed = 42)
  m <- lingam_var(s$data, lags = 1, reg_method = "ols", criterion = NULL, prune = FALSE)

  expected_order <- 1:3
  expected_B0 <- matrix(
    c(0, 0.6270137, 0.0199364,
      0, 0, -0.5048833,
      0, 0, 0),
    nrow = 3, ncol = 3
  )
  expected_B1 <- matrix(
    c(0.3407105, 0.0177596, 0.0553924,
      0.0505192, 0.2804865, -0.075001,
      0.3143521, -0.0262389, 0.4227622),
    nrow = 3, ncol = 3
  )

  expect_equal(m$causal_order, expected_order)
  expect_equal(unname(m$adjacency_matrices[1, , ]), expected_B0, tolerance = 1e-6)
  expect_equal(unname(m$adjacency_matrices[2, , ]), expected_B1, tolerance = 1e-6)
})

test_that("lingam_var snapshot is reproducible across repeated runs", {
  s <- generate_varlingam_sample(n = 500, seed = 42)
  m1 <- lingam_var(s$data, lags = 1, reg_method = "ols", criterion = NULL, prune = FALSE)
  m2 <- lingam_var(s$data, lags = 1, reg_method = "ols", criterion = NULL, prune = FALSE)
  expect_identical(m1$adjacency_matrices, m2$adjacency_matrices)
})
