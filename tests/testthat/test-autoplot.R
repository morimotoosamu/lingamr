test_that("autoplot.LingamResult returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("igraph")

  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  res <- lingam_direct(dat$data, reg_method = "ols")
  pl  <- ggplot2::autoplot(res)

  expect_s3_class(pl, "ggplot")
})

test_that("autoplot.LingamResult works when all edges are filtered out", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("igraph")

  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  res <- lingam_direct(dat$data, reg_method = "ols")
  pl  <- ggplot2::autoplot(res, threshold = 100)  # all edges filtered out

  expect_s3_class(pl, "ggplot")
})

test_that("autoplot.LingamResult respects label_edges = FALSE", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("igraph")

  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  res <- lingam_direct(dat$data, reg_method = "ols")
  pl  <- ggplot2::autoplot(res, label_edges = FALSE)

  expect_s3_class(pl, "ggplot")
})

test_that("autoplot.LingamResult works without column names", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("igraph")

  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  mat <- unname(as.matrix(dat$data))
  res <- lingam_direct(mat, reg_method = "ols")
  pl  <- ggplot2::autoplot(res)

  expect_s3_class(pl, "ggplot")
})
