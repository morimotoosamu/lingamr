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
