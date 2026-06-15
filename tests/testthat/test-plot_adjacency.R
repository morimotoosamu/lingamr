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

  # debug = TRUE outputs the DOT string. Since color names are converted to hex internally,
  # we check the comparison-mode characteristic that each edge row gets a color attribute.
  dot_output <- capture.output(
    plot_adjacency(res$adjacency_matrix,
                   true_B = dat$true_adjacency,
                   debug = TRUE)
  )

  edge_rows <- grep("->", dot_output, value = TRUE)
  expect_true(length(edge_rows) > 0)
  # In comparison mode each edge gets color = '#...' (in normal mode it does not)
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
  # Normal-mode edge rows have label only, with no individual color attribute
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
