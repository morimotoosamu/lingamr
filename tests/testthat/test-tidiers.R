test_that("tidy.LingamResult returns long data.frame of edges", {
  dat <- generate_lingam_sample_6(n = 500, seed = 42)
  res <- lingam_direct(dat$data, reg_method = "ols")
  td  <- tidy(res)

  expect_s3_class(td, "data.frame")
  expect_named(td, c("from", "to", "estimate"))
  expect_type(td$from, "character")
  expect_type(td$to, "character")
  expect_type(td$estimate, "double")
  # 行数 == 非ゼロエッジ数
  expect_equal(nrow(td), sum(abs(res$adjacency_matrix) > 0))
})

test_that("tidy.LingamResult threshold filters edges", {
  dat <- generate_lingam_sample_6(n = 500, seed = 42)
  res <- lingam_direct(dat$data, reg_method = "ols")

  td_all  <- tidy(res, threshold = 0)
  td_high <- tidy(res, threshold = 100)  # 全エッジ除外

  expect_gte(nrow(td_all), nrow(td_high))
  expect_equal(nrow(td_high), 0L)
})

test_that("tidy.LingamResult from/to follow j -> i convention", {
  # 真の構造を直接渡して規則を確認: B["x0","x3"] = 3 は x3 -> x0
  dat  <- generate_lingam_sample_6(n = 100, seed = 1)
  fake <- structure(
    list(adjacency_matrix = dat$true_adjacency, causal_order = 1:6),
    class = "LingamResult"
  )
  td <- tidy(fake)

  edge <- td[td$from == "x3" & td$to == "x0", ]
  expect_equal(nrow(edge), 1L)
  expect_equal(edge$estimate, 3.0)
})

test_that("tidy.LingamResult returns 0-row data.frame for empty graph", {
  fake <- structure(
    list(adjacency_matrix = matrix(0, 3, 3,
           dimnames = list(c("a", "b", "c"), c("a", "b", "c"))),
         causal_order = 1:3),
    class = "LingamResult"
  )
  td <- tidy(fake)

  expect_equal(nrow(td), 0L)
  expect_named(td, c("from", "to", "estimate"))
})

test_that("glance.LingamResult returns one-row summary", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  res <- lingam_direct(dat$data, reg_method = "ols")
  g   <- glance(res)

  expect_s3_class(g, "data.frame")
  expect_equal(nrow(g), 1L)
  expect_true(all(c("n_variables", "n_edges", "causal_order") %in% names(g)))
  expect_equal(g$n_variables, 6L)
})

test_that("tidy.BootstrapResult returns direction counts", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 15L, seed = 42L)
  td  <- tidy(bs)

  expect_s3_class(td, "data.frame")
  expect_true(all(c("from", "to", "count", "proportion") %in% names(td)))
})
