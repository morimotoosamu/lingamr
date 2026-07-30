# Snapshot (golden-value) regression test for lingam_varma().
#
# These pinned numbers are NOT a Python-parity check; they are this R version's
# own deterministic output, frozen to catch *unintended* numerical changes in
# the VARMA-LiNGAM pipeline. The configuration below is fully deterministic
# (uniform-seeded data, Hannan-Rissanen estimation, zero-initialized residual
# filtering, pwling measure, OLS regression, no pruning, no order selection),
# so the result is reproducible across runs.
#
# If the algorithm is changed on purpose, regenerate the expected values:
#   s <- generate_varmalingam_sample(n = 500, seed = 42)
#   m <- lingam_varma(s$data, order = c(1, 1), reg_method = "ols",
#                     criterion = NULL, prune = FALSE)
#   dput(round(unname(m$adjacency_matrices$psis[1, , ]), 7))   # B0
#   dput(round(unname(m$adjacency_matrices$psis[2, , ]), 7))   # psi1
#   dput(round(unname(m$adjacency_matrices$omegas[1, , ]), 7)) # omega1

test_that("lingam_varma output matches the pinned snapshot", {
  s <- generate_varmalingam_sample(n = 500, seed = 42)
  m <- lingam_varma(s$data,
    order = c(1, 1), reg_method = "ols",
    criterion = NULL, prune = FALSE
  )

  expected_order <- 1:3
  expected_B0 <- matrix(
    c(0, 0.6342451, 0.0209408,
      0, 0, -0.5046732,
      0, 0, 0),
    nrow = 3, ncol = 3
  )
  expected_psi1 <- matrix(
    c(0.1528692, -0.1277368, -0.0544168,
      0.0574225, 0.2847738, 0.172352,
      0.1964342, -0.1266163, 0.4224544),
    nrow = 3, ncol = 3
  )
  expected_omega1 <- matrix(
    c(0.3983782, 0.0947546, 0.1275092,
      -0.0058792, 0.1493414, -0.0876858,
      0.0563778, -0.0172027, 0.0301755),
    nrow = 3, ncol = 3
  )

  expect_equal(m$causal_order, expected_order)
  expect_equal(unname(m$adjacency_matrices$psis[1, , ]), expected_B0, tolerance = 1e-6)
  expect_equal(unname(m$adjacency_matrices$psis[2, , ]), expected_psi1, tolerance = 1e-6)
  expect_equal(unname(m$adjacency_matrices$omegas[1, , ]), expected_omega1, tolerance = 1e-6)
})

test_that("lingam_varma snapshot is reproducible across repeated runs", {
  # Direct check that the zero-initialized residual filtering keeps the whole
  # pipeline deterministic (the Python reference is not, because it draws the
  # initial residuals from a standard normal).
  s <- generate_varmalingam_sample(n = 500, seed = 42)
  m1 <- lingam_varma(s$data,
    order = c(1, 1), reg_method = "ols",
    criterion = NULL, prune = FALSE
  )
  m2 <- lingam_varma(s$data,
    order = c(1, 1), reg_method = "ols",
    criterion = NULL, prune = FALSE
  )
  expect_identical(m1$adjacency_matrices, m2$adjacency_matrices)
})
