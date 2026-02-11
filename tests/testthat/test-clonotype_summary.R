test_that("clonotype_summary returns expected structure", {
  vdj <- generate_test_VDJ(n_cells = 30, n_clones = 5, seed = 500)
  result <- assign_clonotype(vdj, num_cores = 1, verbose = FALSE)

  summary_df <- clonotype_summary(result$TCR_data)

  expect_s3_class(summary_df, "data.frame")
  expect_true("n_clones" %in% colnames(summary_df))

  # Should have Match columns for each stringency level
  expected_match_cols <- c(
    "Match.dual_chain_dual_allele",
    "Match.dual_chain_one_partial",
    "Match.dual_chain_both_partial",
    "Match.single_chain_dual_allele",
    "Match.single_chain_single_allele"
  )
  for (col in expected_match_cols) {
    expect_true(col %in% colnames(summary_df),
                info = paste("Missing column:", col))
  }
})

test_that("clonotype_summary errors on missing clone columns", {
  df <- data.frame(x = 1:5)
  expect_error(clonotype_summary(df), "Missing required clone ID columns")
})
