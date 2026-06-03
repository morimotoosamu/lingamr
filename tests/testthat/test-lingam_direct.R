test_that("lingam_direct returns LingamResult with correct structure", {
  dat <- generate_lingam_sample_6(n = 500, seed = 42)
  res <- lingam_direct(dat$data, reg_method = "ols")

  expect_s3_class(res, "LingamResult")
  expect_named(res, c("adjacency_matrix", "causal_order"))
  expect_true(is.matrix(res$adjacency_matrix))
  expect_equal(dim(res$adjacency_matrix), c(6L, 6L))
  expect_equal(length(res$causal_order), 6L)
  # 列名が保持されている
  expect_equal(colnames(res$adjacency_matrix), names(dat$data))
})

test_that("lingam_direct identifies x3 as root (first in causal order)", {
  # x3 は親を持たない外生変数なので因果順序の先頭になるはず
  dat <- generate_lingam_sample_6(n = 2000, seed = 42)
  res <- lingam_direct(dat$data, reg_method = "ols")

  expect_equal(res$causal_order[1], which(names(dat$data) == "x3"))
})

test_that("lingam_direct accepts data.frame and preserves colnames", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  res <- lingam_direct(dat$data, reg_method = "ols")

  expect_equal(colnames(res$adjacency_matrix), names(dat$data))
  expect_equal(rownames(res$adjacency_matrix), names(dat$data))
})

test_that("lingam_direct accepts matrix input", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  mat <- as.matrix(dat$data)
  res <- lingam_direct(mat, reg_method = "ols")

  expect_s3_class(res, "LingamResult")
})

test_that("lingam_direct errors on invalid inputs", {
  dat <- generate_lingam_sample_6(n = 100, seed = 1)

  expect_error(lingam_direct(dat$data, measure = "bad_measure"))
  expect_error(lingam_direct(dat$data, reg_method = "bad_method"))
  expect_error(lingam_direct(dat$data, lambda = "bad_lambda"))
  expect_error(lingam_direct(matrix(1:4, nrow = 2)))   # 非 numeric（整数はOKなので別ケース）
  expect_error(lingam_direct(as.data.frame(matrix(1, nrow = 5, ncol = 1))))  # 1変数
  expect_error(lingam_direct(matrix(numeric(0), nrow = 0, ncol = 3)))        # 0行
})

test_that("print.LingamResult runs without error", {
  dat <- generate_lingam_sample_6(n = 200, seed = 1)
  res <- lingam_direct(dat$data, reg_method = "ols")

  expect_output(print(res), "Direct LiNGAM Result")
  expect_output(print(res), "Causal order")
  expect_output(print(res), "Adjacency matrix")
})

test_that("lingam_direct with prior_knowledge runs without error", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  pk <- make_prior_knowledge(6,
    exogenous_variables = 4,
    labels = names(dat$data)
  )
  res <- lingam_direct(dat$data, prior_knowledge = pk, reg_method = "ols")

  expect_s3_class(res, "LingamResult")
})
