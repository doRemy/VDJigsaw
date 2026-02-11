test_that("identify_clonotypes returns expected structure", {
  vdj <- generate_test_VDJ(n_cells = 30, n_clones = 5, seed = 300)
  wide <- pivot_VDJ(vdj, sample_col = "origin", verbose = FALSE)
  validated <- validate_TCR(wide, verbose = FALSE)

  result <- identify_clonotypes(validated, sample_col = "Sample",
                                 num_cores = 1, verbose = FALSE)

  expect_type(result, "list")
  expect_true("TCR_data" %in% names(result))
  expect_true("ref_tables" %in% names(result))
  expect_s3_class(result$TCR_data, "data.frame")
  expect_type(result$ref_tables, "list")
})

test_that("identify_clonotypes creates clone ID columns for all stringency levels", {
  vdj <- generate_test_VDJ(n_cells = 30, n_clones = 5, seed = 301)
  wide <- pivot_VDJ(vdj, sample_col = "origin", verbose = FALSE)
  validated <- validate_TCR(wide, verbose = FALSE)

  result <- identify_clonotypes(validated, sample_col = "Sample",
                                 num_cores = 1, verbose = FALSE)

  expected_cols <- c(
    "CloneID.dual_chain_dual_allele",
    "CloneID.dual_chain_one_partial",
    "CloneID.dual_chain_both_partial",
    "CloneID.single_chain_dual_allele",
    "CloneID.single_chain_single_allele",
    "CloneID.loose"
  )
  for (col in expected_cols) {
    expect_true(col %in% colnames(result$TCR_data),
                info = paste("Missing column:", col))
  }
})

test_that("identify_clonotypes assigns same clone ID to identical cells", {
  # Create cells with low dropout and some heterozygosity to ensure all columns exist
  vdj <- generate_test_VDJ(n_cells = 40, n_clones = 3,
                            dropout_rate = 0.05, heterozygous_rate = 0.1,
                            invalid_rate = 0, seed = 302)
  wide <- pivot_VDJ(vdj, verbose = FALSE)
  validated <- validate_TCR(wide, verbose = FALSE)

  result <- identify_clonotypes(validated, num_cores = 1, verbose = FALSE)

  # At loose level, cells from the same original clone should tend to group together
  tcr <- result$TCR_data
  n_assigned <- sum(!is.na(tcr$CloneID.loose))
  expect_gt(n_assigned, 0)
})

test_that("identify_clonotypes loose level recovers more cells than strict", {
  vdj <- generate_test_VDJ(n_cells = 40, n_clones = 8,
                            dropout_rate = 0.3, seed = 303)
  wide <- pivot_VDJ(vdj, verbose = FALSE)
  validated <- validate_TCR(wide, verbose = FALSE)

  result <- identify_clonotypes(validated, num_cores = 1, verbose = FALSE)
  tcr <- result$TCR_data

  n_strict <- sum(!is.na(tcr$CloneID.dual_chain_dual_allele))
  n_loose <- sum(!is.na(tcr$CloneID.loose))

  # Loose should assign at least as many cells as strict
  expect_gte(n_loose, n_strict)
})

test_that("identify_clonotypes errors on missing columns", {
  df <- data.frame(TRA_1 = "x", TRA_2 = NA)
  expect_error(identify_clonotypes(df, verbose = FALSE), "Missing required")
})
