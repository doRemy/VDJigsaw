test_that("generate_test_VDJ produces valid output", {
  vdj <- generate_test_VDJ(n_cells = 30, n_clones = 5, seed = 123)

  expect_s3_class(vdj, "data.frame")
  expect_true(nrow(vdj) > 0)

  # Check all required 10x columns exist
  expected_cols <- c(
    "barcode", "is_cell", "contig_id", "high_confidence", "length",
    "chain", "v_gene", "d_gene", "j_gene", "c_gene",
    "full_length", "productive",
    "fwr1", "fwr1_nt", "cdr1", "cdr1_nt",
    "fwr2", "fwr2_nt", "cdr2", "cdr2_nt",
    "fwr3", "fwr3_nt", "cdr3", "cdr3_nt",
    "fwr4", "fwr4_nt",
    "reads", "umis",
    "raw_clonotype_id", "raw_consensus_id", "exact_subclonotype_id",
    "origin"
  )
  for (col in expected_cols) {
    expect_true(col %in% colnames(vdj), info = paste("Missing column:", col))
  }
})

test_that("generate_test_VDJ respects parameters", {
  vdj <- generate_test_VDJ(n_cells = 10, n_clones = 3,
                            samples = c("S1"), seed = 99)

  expect_true(all(vdj$origin == "S1"))
  expect_true(all(vdj$chain %in% c("TRA", "TRB")))
})

test_that("generate_test_VDJ is reproducible with seed", {
  vdj1 <- generate_test_VDJ(n_cells = 20, seed = 42)
  vdj2 <- generate_test_VDJ(n_cells = 20, seed = 42)

  expect_identical(vdj1, vdj2)
})

test_that("generate_test_VDJ includes multiple barcodes per clone", {
  vdj <- generate_test_VDJ(n_cells = 50, n_clones = 5, seed = 42)

  # With 50 cells and 5 clones, at least some clones should have multiple cells
  clones_per_barcode <- table(vdj$raw_clonotype_id)
  expect_true(any(clones_per_barcode > 1))
})
