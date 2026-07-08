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

test_that("autoplot draws Parce / RCD results with NA entries", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("igraph")

  B <- matrix(0, 3, 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  B["b", "a"] <- 1.5
  B["a", "c"] <- NA
  B["c", "a"] <- NA

  parce <- structure(
    list(adjacency_matrix = B, causal_order = list(2L, c(1L, 3L)),
         p_values = matrix(0, 3, 3), independence = "hsic"),
    class = "ParceLingamResult"
  )
  expect_s3_class(ggplot2::autoplot(parce), "ggplot")

  rcd <- structure(
    list(adjacency_matrix = B,
         ancestors_list = list(integer(0), 1L, integer(0))),
    class = "RCDResult"
  )
  expect_s3_class(ggplot2::autoplot(rcd), "ggplot")
})

test_that("autoplot.LiMResult returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("igraph")

  B <- matrix(0, 3, 3, dimnames = list(paste0("x", 1:3), paste0("x", 1:3)))
  B["x2", "x1"] <- 1.2
  lim <- structure(
    list(adjacency_matrix = B, causal_order = 1:3,
         is_continuous = c(TRUE, FALSE, TRUE)),
    class = "LiMResult"
  )

  expect_s3_class(ggplot2::autoplot(lim), "ggplot")
})

test_that("autoplot.MultiGroupLingamResult plots the selected group", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("igraph")

  make_B <- function(coef) {
    B <- matrix(0, 3, 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
    B["b", "a"] <- coef
    B
  }
  fake <- structure(
    list(adjacency_matrices = list(g1 = make_B(1), g2 = make_B(2)),
         causal_order = 1:3),
    class = "MultiGroupLingamResult"
  )

  expect_s3_class(ggplot2::autoplot(fake), "ggplot")
  pl <- ggplot2::autoplot(fake, group = "g2")
  expect_s3_class(pl, "ggplot")
  expect_identical(pl$labels$subtitle, "g2")
})
