test_that("assign_clonotype runs full pipeline end-to-end", {
  vdj <- generate_test_VDJ(n_cells = 30, n_clones = 5, seed = 400)

  result <- assign_clonotype(vdj, sample_col = "origin",
                              num_cores = 1, verbose = FALSE)

  expect_type(result, "list")
  expect_true("TCR_data" %in% names(result))
  expect_true("ref_tables" %in% names(result))
  expect_s3_class(result$TCR_data, "data.frame")
})

test_that("assign_clonotype output has CloneID column", {
  vdj <- generate_test_VDJ(n_cells = 20, n_clones = 5, seed = 401)

  result <- assign_clonotype(vdj, sample_col = "origin",
                              num_cores = 1, verbose = FALSE)

  expect_true("CloneID" %in% colnames(result$TCR_data))
  expect_true("CloneID.loose" %in% colnames(result$TCR_data))
})

test_that("assign_clonotype assigns at least some cells", {
  vdj <- generate_test_VDJ(n_cells = 30, n_clones = 5, seed = 402)

  result <- assign_clonotype(vdj, sample_col = "origin",
                              num_cores = 1, verbose = FALSE)

  n_assigned <- sum(!is.na(result$TCR_data$CloneID))
  expect_gt(n_assigned, 0)
})

test_that("assign_clonotype works without sample_col", {
  vdj <- generate_test_VDJ(n_cells = 20, n_clones = 5, seed = 403)

  result <- assign_clonotype(vdj, sample_col = NULL,
                              num_cores = 1, verbose = FALSE)

  expect_type(result, "list")
  expect_true(nrow(result$TCR_data) > 0)
})

test_that("assign_clonotype ref_tables have expected stringency levels", {
  vdj <- generate_test_VDJ(n_cells = 20, n_clones = 5, seed = 404)

  result <- assign_clonotype(vdj, num_cores = 1, verbose = FALSE)

  expected_levels <- c(
    "dual_chain_dual_allele", "dual_chain_one_partial",
    "dual_chain_both_partial", "single_chain_dual_allele",
    "single_chain_single_allele", "loose"
  )
  for (level in expected_levels) {
    expect_true(level %in% names(result$ref_tables),
                info = paste("Missing ref_table level:", level))
  }
})
