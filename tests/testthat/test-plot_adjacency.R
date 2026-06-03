test_that("plot_adjacency errors without DiagrammeR", {
  skip_if(requireNamespace("DiagrammeR", quietly = TRUE), "DiagrammeR is installed")
  dat <- generate_lingam_sample_6(n = 100, seed = 1)
  expect_error(plot_adjacency(dat$true_adjacency))
})

test_that("plot_adjacency returns grViz object in normal mode", {
  skip_if_not_installed("DiagrammeR")
  dat <- generate_lingam_sample_6(n = 100, seed = 1)
  g <- plot_adjacency(dat$true_adjacency)
  expect_s3_class(g, "grViz")
})

test_that("plot_adjacency returns grViz object in comparison mode", {
  skip_if_not_installed("DiagrammeR")
  dat <- generate_lingam_sample_6(n = 500, seed = 42)
  res <- lingam_direct(dat$data, reg_method = "ols")
  g <- plot_adjacency(res$adjacency_matrix, true_B = dat$true_adjacency)
  expect_s3_class(g, "grViz")
})

test_that("plot_adjacency comparison mode adds per-edge color attributes", {
  skip_if_not_installed("DiagrammeR")
  dat <- generate_lingam_sample_6(n = 500, seed = 42)
  res <- lingam_direct(dat$data, reg_method = "ols")

  # debug = TRUE で DOT 文字列を出力させる。色名は内部で hex に変換されるため、
  # 比較モードの特徴である「各エッジ行に color 属性が付くこと」を確認する。
  dot_output <- capture.output(
    plot_adjacency(res$adjacency_matrix,
                   true_B = dat$true_adjacency,
                   debug = TRUE)
  )

  edge_rows <- grep("->", dot_output, value = TRUE)
  expect_true(length(edge_rows) > 0)
  # 比較モードではエッジごとに color = '#...' が付与される（通常モードでは付かない）
  expect_true(any(grepl("color = '#", edge_rows)))
})

test_that("plot_adjacency normal mode does NOT add per-edge color attributes", {
  skip_if_not_installed("DiagrammeR")
  dat <- generate_lingam_sample_6(n = 500, seed = 42)
  res <- lingam_direct(dat$data, reg_method = "ols")

  dot_output <- capture.output(
    plot_adjacency(res$adjacency_matrix, debug = TRUE)
  )

  edge_rows <- grep("->", dot_output, value = TRUE)
  expect_true(length(edge_rows) > 0)
  # 通常モードのエッジ行は label のみで、個別 color 属性は付かない
  expect_false(any(grepl("color =", edge_rows)))
})

test_that("plot_adjacency errors when true_B has wrong dimensions", {
  skip_if_not_installed("DiagrammeR")
  dat <- generate_lingam_sample_6(n = 100, seed = 1)
  bad_true_B <- matrix(0, nrow = 3, ncol = 3)
  expect_error(
    plot_adjacency(dat$true_adjacency, true_B = bad_true_B),
    "same dimensions"
  )
})

test_that("plot_adjacency returns invisible NULL when no edges above threshold", {
  skip_if_not_installed("DiagrammeR")
  B <- matrix(0.01, nrow = 3, ncol = 3)
  diag(B) <- 0
  expect_invisible(plot_adjacency(B, threshold = 1.0))
})

test_that("plot_adjacency rankdir validation works", {
  skip_if_not_installed("DiagrammeR")
  dat <- generate_lingam_sample_6(n = 100, seed = 1)
  expect_error(plot_adjacency(dat$true_adjacency, rankdir = "XX"))
})
