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
  # 既知の係数
  expect_equal(B["x0", "x3"],  3.0)
  expect_equal(B["x2", "x3"],  6.0)
  expect_equal(B["x1", "x0"],  3.0)
  expect_equal(B["x1", "x2"],  2.0)
  expect_equal(B["x5", "x0"],  4.0)
  expect_equal(B["x4", "x0"],  8.0)
  expect_equal(B["x4", "x2"], -1.0)
  # 真の構造にない辺はゼロ
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
