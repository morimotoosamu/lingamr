# ── input validation ──────────────────────────────────────────────────────────

test_that("lingam_varma validates inputs", {
  X <- varmas_1000_s42()$data

  expect_error(lingam_varma("abc"), "numeric matrix")
  expect_error(lingam_varma(X[, 1, drop = FALSE]), "at least 2 variables")
  expect_error(lingam_varma(X[1:2, ]), "at least 3 observations")
  Xna <- X
  Xna[1, 1] <- NA
  expect_error(lingam_varma(Xna), "missing values")

  expect_error(lingam_varma(X, order = c(1, 1, 1)), "order must be")
  expect_error(lingam_varma(X, order = c(-1, 1)), "order must be")
  expect_error(lingam_varma(X, order = c(0, 0)), "order must be")
  expect_error(lingam_varma(X, order = c(1, 1), criterion = "fpe"), "arg")
  expect_error(lingam_varma(X, order = c(1, 1), prune = NA), "prune must be")
})

test_that("lingam_varma validates ar_coefs / ma_coefs", {
  X <- varmas_1000_s42()$data
  phi <- array(0.3 * diag(3), c(1, 3, 3))
  theta <- array(0.2 * diag(3), c(1, 3, 3))

  expect_error(lingam_varma(X, ar_coefs = phi), "supplied together")
  expect_error(lingam_varma(X, ma_coefs = theta), "supplied together")
  expect_error(
    lingam_varma(X, ar_coefs = matrix(0, 3, 3), ma_coefs = theta),
    "numeric array"
  )
  expect_error(
    lingam_varma(X, ar_coefs = array(0, c(1, 2, 2)), ma_coefs = theta),
    "shape"
  )
})

# ── structure ─────────────────────────────────────────────────────────────────

test_that("lingam_varma returns the documented structure", {
  m <- fit_varma_1000()

  expect_s3_class(m, "VARMALiNGAMResult")
  expect_named(m, c(
    "adjacency_matrices", "causal_order", "residuals", "order",
    "ar_coefs", "ma_coefs", "const"
  ))
  expect_equal(m$order, c(1L, 1L))

  psis <- m$adjacency_matrices$psis
  omegas <- m$adjacency_matrices$omegas
  expect_equal(dim(psis), c(2L, 3L, 3L))
  expect_equal(dim(omegas), c(1L, 3L, 3L))
  expect_equal(dimnames(psis)[[1]], c("lag0", "lag1"))
  expect_equal(dimnames(psis)[[2]], c("x0", "x1", "x2"))
  expect_equal(dimnames(omegas)[[1]], "ma1")

  # residuals drop the first max(p, q) rows
  expect_equal(dim(m$residuals), c(999L, 3L))
  expect_equal(colnames(m$residuals), c("x0", "x1", "x2"))

  # B0 is acyclic with a zero diagonal: permuted by the causal order it is
  # strictly lower triangular
  B0 <- psis[1, , ]
  expect_true(all(diag(B0) == 0))
  perm <- B0[m$causal_order, m$causal_order]
  expect_true(all(perm[upper.tri(perm, diag = TRUE)] == 0))

  expect_equal(dim(m$ar_coefs), c(1L, 3L, 3L))
  expect_equal(dim(m$ma_coefs), c(1L, 3L, 3L))
  expect_length(m$const, 3L)
})

# ── coefficient recovery ──────────────────────────────────────────────────────

test_that("lingam_varma recovers the true structure on n = 2000", {
  s <- varmas_2000_s42()
  m <- fit_varma_2000()

  expect_equal(m$causal_order, 1:3)
  expect_lt(max(abs(m$adjacency_matrices$psis[1, , ] - s$true_B0)), 0.1)
  expect_lt(max(abs(m$adjacency_matrices$psis[2, , ] - s$true_psi1)), 0.1)
  expect_lt(max(abs(m$adjacency_matrices$omegas[1, , ] - s$true_omega1)), 0.15)
})

test_that("unpruned omegas equal the exact similarity transform", {
  # Permanent guard for the upstream bug fix: the Python reference computes
  # (I - B0) Theta only (three-argument np.dot treats the inverse as an output
  # buffer); this implementation applies the full transform
  # omega = (I - B0) Theta (I - B0)^{-1}.
  m <- fit_varma_1000()
  B0 <- m$adjacency_matrices$psis[1, , ]
  ib0 <- diag(3) - B0
  expected <- ib0 %*% m$ma_coefs[1, , ] %*% solve(ib0)
  expect_equal(unname(m$adjacency_matrices$omegas[1, , ]), unname(expected),
    tolerance = 1e-12
  )
  # psi1 likewise equals (I - B0) Phi1
  expect_equal(
    unname(m$adjacency_matrices$psis[2, , ]),
    unname(ib0 %*% m$ar_coefs[1, , ]),
    tolerance = 1e-12
  )
})

# ── degenerate orders ─────────────────────────────────────────────────────────

test_that("lingam_varma handles pure AR and pure MA orders", {
  X <- varmas_1000_s42()$data

  m10 <- lingam_varma(X,
    order = c(1, 0), criterion = NULL,
    reg_method = "ols", prune = FALSE
  )
  expect_equal(m10$order, c(1L, 0L))
  expect_equal(dim(m10$adjacency_matrices$psis), c(2L, 3L, 3L))
  expect_equal(dim(m10$adjacency_matrices$omegas), c(0L, 3L, 3L))
  expect_equal(dim(m10$residuals), c(999L, 3L))

  m01 <- lingam_varma(X,
    order = c(0, 1), criterion = NULL,
    reg_method = "ols", prune = FALSE
  )
  expect_equal(m01$order, c(0L, 1L))
  expect_equal(dim(m01$adjacency_matrices$psis), c(1L, 3L, 3L))
  expect_equal(dim(m01$adjacency_matrices$omegas), c(1L, 3L, 3L))
})

# ── order selection ───────────────────────────────────────────────────────────

test_that("order selection returns a valid order and is reflected in the result", {
  X <- varmas_2000_s42()$data
  m <- lingam_varma(X,
    order = c(2, 2), criterion = "bic",
    reg_method = "ols", prune = FALSE
  )
  # A VARMA(1,1) DGP is well approximated by a slightly longer pure VAR, so
  # the selected order is not pinned; require a valid, non-degenerate choice.
  expect_true(all(m$order >= 0L) && all(m$order <= 2L) && sum(m$order) >= 1L)
  expect_equal(dim(m$adjacency_matrices$psis), c(1L + m$order[1], 3L, 3L))
  expect_equal(dim(m$adjacency_matrices$omegas), c(m$order[2], 3L, 3L))
})

test_that("Hannan-Rissanen guards against too small samples", {
  # far too small: the stage-1 long VAR is infeasible even after shrinking h
  X12 <- as.matrix(varmas_1000_s42()$data[1:12, ])
  expect_error(
    select_varma_order(X12, max_p = 2L, max_q = 2L, criterion = "bic"),
    "Not enough observations"
  )

  # small but feasible: h is reduced with a warning
  X30 <- as.matrix(varmas_1000_s42()$data[1:30, ])
  expect_warning(hr <- fit_varma_hr(X30, 2L, 2L), "reduced from")
  expect_lt(hr$h, 10L)
})

# ── known coefficients ────────────────────────────────────────────────────────

test_that("lingam_varma skips estimation when ar_coefs and ma_coefs are given", {
  s <- varmas_1000_s42()
  m <- lingam_varma(s$data,
    ar_coefs = array(s$true_phi1, c(1, 3, 3)),
    ma_coefs = array(s$true_theta1, c(1, 3, 3)),
    criterion = NULL, reg_method = "ols", prune = FALSE
  )
  expect_equal(m$ar_coefs, array(s$true_phi1, c(1, 3, 3)))
  expect_equal(m$ma_coefs, array(s$true_theta1, c(1, 3, 3)))
  expect_equal(m$const, numeric(3))
  expect_equal(m$causal_order, 1:3)
  expect_lt(max(abs(m$adjacency_matrices$psis[1, , ] - s$true_B0)), 0.1)
})

# ── pruning ───────────────────────────────────────────────────────────────────

test_that("pruning shrinks weak edges to zero", {
  skip_if_not_installed("glmnet")
  s <- varmas_1000_s42()
  m <- lingam_varma(s$data,
    order = c(1, 1), criterion = NULL,
    reg_method = "ols", prune = TRUE
  )
  psis <- m$adjacency_matrices$psis
  # true zero entries of B0 should be sparse after pruning
  expect_gt(sum(psis[1, , ] == 0), 6L) # 7 true zeros in the 3x3 B0
  # the strong true edges survive
  expect_gt(abs(psis[1, 2, 1]), 0.3)
  expect_gt(abs(psis[1, 3, 2]), 0.2)
})

# ── invertibility warning ─────────────────────────────────────────────────────

test_that("non-invertible ma_coefs trigger a warning", {
  # A mildly explosive MA on a short series: the warning fires while the
  # residual recursion stays finite enough for the rest of the pipeline.
  s <- varmas_1000_s42()
  expect_warning(
    lingam_varma(s$data[1:100, ],
      ar_coefs = array(s$true_phi1, c(1, 3, 3)),
      ma_coefs = array(1.05 * diag(3), c(1, 3, 3)),
      criterion = NULL, reg_method = "ols", prune = FALSE
    ),
    "not invertible"
  )
})

# ── print ─────────────────────────────────────────────────────────────────────

test_that("print.VARMALiNGAMResult shows the key components", {
  m <- fit_varma_1000()
  out <- capture.output(print(m))
  expect_true(any(grepl("VARMA-LiNGAM Result", out)))
  expect_true(any(grepl("Order \\(p, q\\) : \\(1, 1\\)", out)))
  expect_true(any(grepl("x0 -> x1 -> x2", out)))
  expect_true(any(grepl("Instantaneous adjacency matrix B0", out)))
  expect_true(any(grepl("psi1", out)))
  expect_true(any(grepl("omega1", out)))
})
