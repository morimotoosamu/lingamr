test_that("generate_lingam_sample_6 returns correct structure", {
  out <- generate_lingam_sample_6(n = 200, seed = 1)

  expect_named(out, c("data", "true_adjacency"))
  expect_s3_class(out$data, "data.frame")
  expect_equal(dim(out$data), c(200L, 6L))
  expect_equal(names(out$data), c("x0", "x1", "x2", "x3", "x4", "x5"))
  expect_false(anyNA(out$data))

  B <- out$true_adjacency
  expect_true(is.matrix(B))
  expect_equal(dim(B), c(6L, 6L))
  # known coefficients
  expect_equal(B["x0", "x3"],  3.0)
  expect_equal(B["x2", "x3"],  6.0)
  expect_equal(B["x1", "x0"],  3.0)
  expect_equal(B["x1", "x2"],  2.0)
  expect_equal(B["x5", "x0"],  4.0)
  expect_equal(B["x4", "x0"],  8.0)
  expect_equal(B["x4", "x2"], -1.0)
  # edges not in the true structure are zero
  expect_equal(B["x3", "x0"],  0.0)
})

test_that("generate_lingam_sample_6 seed is reproducible", {
  a <- generate_lingam_sample_6(n = 100, seed = 7)
  b <- generate_lingam_sample_6(n = 100, seed = 7)
  expect_equal(a$data, b$data)
})

test_that("generate_lingam_sample_10 returns correct dimensions", {
  out <- generate_lingam_sample_10(n = 300, seed = 1)

  expect_equal(dim(out$data), c(300L, 10L))
  expect_equal(dim(out$true_adjacency), c(10L, 10L))
  expect_false(anyNA(out$data))
})

test_that("generate_lingam_paradox_data returns correct structure and true chain", {
  out <- generate_lingam_paradox_data(n = 500, seed = 1)

  expect_named(out, c("data", "true_adjacency"))
  expect_s3_class(out$data, "data.frame")
  expect_equal(dim(out$data), c(500L, 4L))
  expect_equal(names(out$data), c("x0", "x1", "x2", "x3"))
  expect_false(anyNA(out$data))

  B <- out$true_adjacency
  expect_true(is.matrix(B))
  expect_equal(dim(B), c(4L, 4L))
  # true causal chain x0 -> x1 -> x2 -> x3, coefficient 0.8 each
  expect_equal(B["x1", "x0"], 0.8)
  expect_equal(B["x2", "x1"], 0.8)
  expect_equal(B["x3", "x2"], 0.8)
  # no other edges
  expect_equal(sum(B != 0), 3L)
})

test_that("generate_lingam_paradox_data seed is reproducible", {
  a <- generate_lingam_paradox_data(n = 200, seed = 3)
  b <- generate_lingam_paradox_data(n = 200, seed = 3)
  expect_equal(a$data, b$data)
})

test_that("generate_lingam_hard_sample returns correct structure", {
  out <- generate_lingam_hard_sample(n = 200, seed = 1)

  expect_named(out, c("data", "true_adjacency"))
  expect_true(ncol(out$data) >= 2L)
  expect_false(anyNA(out$data))
})

test_that("generate_lingam_sample_6 noise_dist variants work", {
  for (dist in c("uniform", "gaussian", "lognormal", "exponential")) {
    out <- generate_lingam_sample_6(n = 50, seed = 1, noise_dist = dist)
    expect_false(anyNA(out$data), label = paste("NA in noise_dist =", dist))
  }
})


# ── generate_lingam_large_sample ──────────────────────────────────────────────

test_that("generate_lingam_large_sample returns correct structure", {
  out <- generate_lingam_large_sample(p = 10, n = 100, seed = 1)

  expect_named(out, c("data", "true_adjacency", "true_causal_order"))
  expect_s3_class(out$data, "data.frame")
  expect_equal(dim(out$data), c(100L, 10L))
  expect_equal(names(out$data), paste0("x", 0:9))
  expect_false(anyNA(out$data))

  B <- out$true_adjacency
  expect_true(is.matrix(B))
  expect_equal(dim(B), c(10L, 10L))
  expect_equal(rownames(B), paste0("x", 0:9))
  expect_equal(colnames(B), paste0("x", 0:9))
})

test_that("generate_lingam_large_sample true_causal_order is 0:(p-1)", {
  out <- generate_lingam_large_sample(p = 8, n = 50, seed = 1)
  expect_equal(out$true_causal_order, 0:7)
})

test_that("generate_lingam_large_sample adjacency matrix is strictly lower-triangular (valid DAG)", {
  out <- generate_lingam_large_sample(p = 15, n = 100, seed = 7)
  B   <- out$true_adjacency
  # upper triangle (including diagonal) is all zero => variables ordered by causal order
  expect_true(all(B[upper.tri(B, diag = TRUE)] == 0))
  # at least one edge exists
  expect_gt(sum(B != 0), 0L)
})

test_that("generate_lingam_large_sample seed is reproducible", {
  a <- generate_lingam_large_sample(p = 8, n = 100, seed = 99)
  b <- generate_lingam_large_sample(p = 8, n = 100, seed = 99)
  expect_equal(a$data,           b$data)
  expect_equal(a$true_adjacency, b$true_adjacency)
})

test_that("generate_lingam_large_sample different seeds give different results", {
  a <- generate_lingam_large_sample(p = 8, n = 100, seed = 1)
  b <- generate_lingam_large_sample(p = 8, n = 100, seed = 2)
  expect_false(identical(a$data, b$data))
})

test_that("generate_lingam_large_sample respects p parameter", {
  for (p_val in c(5L, 15L, 25L)) {
    out <- generate_lingam_large_sample(p = p_val, n = 50, seed = 1)
    expect_equal(ncol(out$data), p_val, label = paste("ncol for p =", p_val))
    expect_equal(nrow(out$data), 50L,   label = paste("nrow for p =", p_val))
    expect_equal(dim(out$true_adjacency), c(p_val, p_val),
                 label = paste("adj dim for p =", p_val))
  }
})

test_that("generate_lingam_large_sample coefficients respect coef range", {
  out   <- generate_lingam_large_sample(p = 10, n = 50,
                                        coef_min = 1.0, coef_max = 2.0, seed = 1)
  edges <- out$true_adjacency[out$true_adjacency != 0]
  expect_gt(length(edges), 0L)                    # at least one edge
  expect_true(all(abs(edges) >= 1.0))
  expect_true(all(abs(edges) <= 2.0))
})

test_that("generate_lingam_large_sample noise_dist variants work", {
  for (dist in c("uniform", "lognormal", "exponential", "t3")) {
    out <- generate_lingam_large_sample(p = 5, n = 50, seed = 1, noise_dist = dist)
    expect_false(anyNA(out$data), label = paste("NA in noise_dist =", dist))
  }
})

test_that("generate_lingam_large_sample input validation", {
  expect_error(generate_lingam_large_sample(p = 1),
               "p must")
  expect_error(generate_lingam_large_sample(n = 1),
               "n must")
  expect_error(generate_lingam_large_sample(max_parents = 0),
               "max_parents")
  expect_error(generate_lingam_large_sample(coef_min = -1),
               "coef_min")
  expect_error(generate_lingam_large_sample(coef_min = 1.5, coef_max = 1.0),
               "coef_min")
  expect_error(generate_lingam_large_sample(noise_dist = "invalid"),
               "noise_dist")
})


# ── generate_lim_sample (Poisson mode) ────────────────────────────────────────

test_that("generate_lim_sample generates Poisson counts when is_poisson = TRUE", {
  dat <- generate_lim_sample(n = 500, seed = 1, is_poisson = TRUE)

  expect_named(dat, c("data", "adjacency_matrix", "is_continuous", "is_poisson"))
  x2 <- dat$data$x2
  expect_true(all(x2 >= 0))
  expect_true(all(x2 == round(x2)))
  expect_gt(max(x2), 1) # actual counts, not binary
  expect_true(dat$is_poisson)
  expect_equal(dat$is_continuous, c(TRUE, FALSE, TRUE))
})

test_that("generate_lim_sample defaults to binary and records is_poisson = FALSE", {
  dat <- generate_lim_sample(n = 200, seed = 1)

  expect_false(dat$is_poisson)
  expect_true(all(dat$data$x2 %in% c(0, 1)))
})

test_that("generate_lim_sample Poisson mode is reproducible and validated", {
  a <- generate_lim_sample(n = 200, seed = 7, is_poisson = TRUE)
  b <- generate_lim_sample(n = 200, seed = 7, is_poisson = TRUE)
  expect_equal(a, b)

  expect_error(generate_lim_sample(n = 100, is_poisson = "yes"), "is_poisson")
  expect_error(generate_lim_sample(n = 100, is_poisson = NA), "is_poisson")
})


# ── generate_resit_sample ─────────────────────────────────────────────────────

test_that("generate_resit_sample returns the documented structure", {
  dat <- generate_resit_sample(n = 100, seed = 1)

  expect_named(dat, c("data", "adjacency_matrix", "causal_order"))
  expect_s3_class(dat$data, "data.frame")
  expect_equal(dim(dat$data), c(100L, 4L))
  expect_equal(names(dat$data), paste0("x", 0:3))

  B <- dat$adjacency_matrix
  expect_equal(dim(B), c(4L, 4L))
  expect_true(all(B %in% c(0, 1)))
  # true edges, m[to, from] convention
  expect_equal(B["x1", "x0"], 1)
  expect_equal(B["x2", "x0"], 1)
  expect_equal(B["x2", "x1"], 1)
  expect_equal(B["x3", "x2"], 1)
  expect_equal(sum(B), 4)

  expect_equal(dat$causal_order, 1:4)
})

test_that("generate_resit_sample is reproducible with a seed and validates inputs", {
  a <- generate_resit_sample(n = 150, seed = 7)
  b <- generate_resit_sample(n = 150, seed = 7)
  expect_equal(a, b)

  expect_error(generate_resit_sample(n = 1), "n must be")
  expect_error(generate_resit_sample(n = 100, seed = "x"), "seed")
})


# ── generate_camuv_sample ─────────────────────────────────────────────────────

test_that("generate_camuv_sample returns the documented structure", {
  dat <- generate_camuv_sample(n = 100, seed = 1)

  expect_named(dat, c("data", "adjacency_matrix", "confounded_pairs"))
  expect_s3_class(dat$data, "data.frame")
  expect_equal(dim(dat$data), c(100L, 6L))
  expect_equal(names(dat$data), paste0("x", 0:5))

  B <- dat$adjacency_matrix
  expect_equal(dim(B), c(6L, 6L))
  expect_equal(dimnames(B), list(names(dat$data), names(dat$data)))
  expect_true(all(B %in% c(0, 1) | is.na(B)))
  expect_equal(B["x1", "x0"], 1)
  expect_equal(B["x3", "x0"], 1)
  expect_equal(B["x4", "x2"], 1)
  expect_equal(sum(B > 0, na.rm = TRUE), 3)

  # UCP pair {x2, x5} and UBP pair {x3, x4} are NA in both directions
  expect_true(is.na(B["x2", "x5"]) && is.na(B["x5", "x2"]))
  expect_true(is.na(B["x3", "x4"]) && is.na(B["x4", "x3"]))
  expect_equal(sum(is.na(B)), 4)

  expect_equal(unname(dat$confounded_pairs), rbind(c(3L, 6L), c(4L, 5L)))
})

test_that("generate_camuv_sample is reproducible with a seed and validates inputs", {
  a <- generate_camuv_sample(n = 150, seed = 7)
  b <- generate_camuv_sample(n = 150, seed = 7)
  expect_equal(a, b)

  expect_error(generate_camuv_sample(n = 1), "n must be")
  expect_error(generate_camuv_sample(n = 100, seed = "x"), "seed")
})


# ── generate_varmalingam_sample ───────────────────────────────────────────────

test_that("generate_varmalingam_sample returns the documented structure", {
  out <- generate_varmalingam_sample(n = 200, seed = 1)

  expect_named(out, c(
    "data", "true_B0", "true_phi1", "true_theta1", "true_psi1", "true_omega1"
  ))
  expect_s3_class(out$data, "data.frame")
  expect_equal(dim(out$data), c(200L, 3L))
  expect_equal(names(out$data), c("x0", "x1", "x2"))
  expect_false(anyNA(out$data))

  # known instantaneous structure x0 -> x1 -> x2
  expect_equal(out$true_B0[2, 1], 0.6)
  expect_equal(out$true_B0[3, 2], -0.5)
  expect_equal(sum(out$true_B0 != 0), 2L)
})

test_that("generate_varmalingam_sample true matrices are internally consistent", {
  out <- generate_varmalingam_sample(n = 50, seed = 1)
  I3 <- diag(3)

  expect_equal(out$true_psi1, (I3 - out$true_B0) %*% out$true_phi1)
  expect_equal(
    out$true_omega1,
    (I3 - out$true_B0) %*% out$true_theta1 %*% solve(I3 - out$true_B0)
  )

  # DGP is stationary (AR eigenvalues inside unit circle) and invertible
  # (MA eigenvalues inside unit circle)
  expect_lt(max(Mod(eigen(out$true_phi1)$values)), 1)
  expect_lt(max(Mod(eigen(out$true_theta1)$values)), 1)
})

test_that("generate_varmalingam_sample seed is reproducible", {
  a <- generate_varmalingam_sample(n = 100, seed = 7)
  b <- generate_varmalingam_sample(n = 100, seed = 7)
  expect_equal(a, b)
  expect_false(identical(
    generate_varmalingam_sample(n = 100, seed = 1)$data,
    generate_varmalingam_sample(n = 100, seed = 2)$data
  ))
})
