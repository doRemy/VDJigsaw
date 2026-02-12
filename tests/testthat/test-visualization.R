# Helper: run the full pipeline to get test data for visualization
.get_test_result <- function(seed = 500) {
  vdj <- generate_test_VDJ(n_cells = 40, n_clones = 8, seed = seed)
  assign_clonotype(vdj, sample_col = "origin", num_cores = 1, verbose = FALSE)
}

# ==============================================================================
# plot_stringency_summary
# ==============================================================================

test_that("plot_stringency_summary returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  result <- .get_test_result()

  p <- plot_stringency_summary(result$TCR_data)
  expect_s3_class(p, "ggplot")
})

test_that("plot_stringency_summary works with sample_filter", {
  skip_if_not_installed("ggplot2")
  result <- .get_test_result()

  p <- plot_stringency_summary(result$TCR_data, sample_filter = "SampleA")
  expect_s3_class(p, "ggplot")
})

test_that("plot_stringency_summary errors on missing columns", {
  skip_if_not_installed("ggplot2")

  bad_data <- data.frame(x = 1:5)
  expect_error(plot_stringency_summary(bad_data), "Missing required clone ID columns")
})

test_that("plot_stringency_summary errors on invalid sample_filter", {
  skip_if_not_installed("ggplot2")
  result <- .get_test_result()

  expect_error(
    plot_stringency_summary(result$TCR_data, sample_filter = "NonExistent"),
    "No data remaining"
  )
})

test_that("plot_stringency_summary per_sample returns list when combined = FALSE", {
  skip_if_not_installed("ggplot2")
  result <- .get_test_result()

  plots <- plot_stringency_summary(result$TCR_data, per_sample = TRUE, combined = FALSE)
  expect_type(plots, "list")
  expect_true(length(plots) > 0)
  for (p in plots) {
    expect_s3_class(p, "ggplot")
  }
})

test_that("plot_stringency_summary per_sample combined returns patchwork", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  result <- .get_test_result()

  p <- plot_stringency_summary(result$TCR_data, per_sample = TRUE, combined = TRUE)
  expect_s3_class(p, "patchwork")
})

# ==============================================================================
# plot_clone_composition
# ==============================================================================

test_that("plot_clone_composition returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  result <- .get_test_result()

  ref <- result$ref_tables$dual_chain_both_partial
  p <- plot_clone_composition(ref)
  expect_s3_class(p, "ggplot")
})

test_that("plot_clone_composition respects top_n", {
  skip_if_not_installed("ggplot2")
  result <- .get_test_result()

  ref <- result$ref_tables$single_chain_single_allele
  p <- plot_clone_composition(ref, top_n = 3)
  expect_s3_class(p, "ggplot")
})

test_that("plot_clone_composition errors on missing columns", {
  skip_if_not_installed("ggplot2")

  bad_data <- data.frame(x = 1:5)
  expect_error(plot_clone_composition(bad_data), "Missing required columns")
})

test_that("plot_clone_composition per_sample returns list when combined = FALSE", {
  skip_if_not_installed("ggplot2")
  result <- .get_test_result()

  ref <- result$ref_tables$dual_chain_both_partial
  # ref_tables may not have Sample column; check before testing
  if ("Sample" %in% colnames(ref)) {
    plots <- plot_clone_composition(ref, per_sample = TRUE, combined = FALSE)
    expect_type(plots, "list")
    expect_true(length(plots) > 0)
    for (p in plots) {
      expect_s3_class(p, "ggplot")
    }
  } else {
    expect_error(
      plot_clone_composition(ref, per_sample = TRUE),
      "Sample column not found"
    )
  }
})

# ==============================================================================
# plot_clonotype_flow
# ==============================================================================

test_that("plot_clonotype_flow returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggalluvial")
  result <- .get_test_result()

  p <- plot_clonotype_flow(result$TCR_data)
  expect_s3_class(p, "ggplot")
})

test_that("plot_clonotype_flow works with sample_filter", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggalluvial")
  result <- .get_test_result()

  p <- plot_clonotype_flow(result$TCR_data, sample_filter = "SampleA")
  expect_s3_class(p, "ggplot")
})

test_that("plot_clonotype_flow works with show_unassigned", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggalluvial")
  result <- .get_test_result()

  p <- plot_clonotype_flow(result$TCR_data, show_unassigned = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_clonotype_flow respects top_n", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggalluvial")
  result <- .get_test_result()

  p <- plot_clonotype_flow(result$TCR_data, top_n = 3)
  expect_s3_class(p, "ggplot")
})

test_that("plot_clonotype_flow errors on missing columns", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggalluvial")

  bad_data <- data.frame(x = 1:5)
  expect_error(plot_clonotype_flow(bad_data), "Missing required clone ID columns")
})

test_that("plot_clonotype_flow errors on invalid sample_filter", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggalluvial")
  result <- .get_test_result()

  expect_error(
    plot_clonotype_flow(result$TCR_data, sample_filter = "NonExistent"),
    "No data remaining"
  )
})

test_that("plot_clonotype_flow per_sample returns list when combined = FALSE", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggalluvial")
  result <- .get_test_result()

  plots <- plot_clonotype_flow(result$TCR_data, per_sample = TRUE, combined = FALSE)
  expect_type(plots, "list")
  expect_true(length(plots) > 0)
  for (p in plots) {
    expect_s3_class(p, "ggplot")
  }
})

test_that("plot_clonotype_flow per_sample combined returns patchwork", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggalluvial")
  skip_if_not_installed("patchwork")
  result <- .get_test_result()

  p <- plot_clonotype_flow(result$TCR_data, per_sample = TRUE, combined = TRUE)
  expect_s3_class(p, "patchwork")
})
