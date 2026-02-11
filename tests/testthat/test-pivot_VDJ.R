test_that("pivot_VDJ produces one row per barcode", {
  vdj <- generate_test_VDJ(n_cells = 30, n_clones = 5, seed = 100)
  wide <- pivot_VDJ(vdj, sample_col = "origin", verbose = FALSE)

  expect_s3_class(wide, "data.frame")
  expect_false(any(duplicated(wide$barcode)))
  expect_equal(nrow(wide), length(unique(vdj$barcode)))
})

test_that("pivot_VDJ creates TRA/TRB chain columns", {
  vdj <- generate_test_VDJ(n_cells = 20, n_clones = 5, seed = 101)
  wide <- pivot_VDJ(vdj, verbose = FALSE)

  expect_true(all(c("TRA_1", "TRA_2", "TRB_1", "TRB_2") %in% colnames(wide)))
  expect_true("Sample" %in% colnames(wide))
  expect_true("barcode" %in% colnames(wide))
})

test_that("pivot_VDJ uses default Sample when sample_col is NULL", {
  vdj <- generate_test_VDJ(n_cells = 10, n_clones = 3, seed = 102)
  wide <- pivot_VDJ(vdj, sample_col = NULL, verbose = FALSE)

  expect_true(all(wide$Sample == "Sample"))
})

test_that("pivot_VDJ errors on invalid sample_col", {
  vdj <- generate_test_VDJ(n_cells = 10, n_clones = 3, seed = 103)
  expect_error(pivot_VDJ(vdj, sample_col = "nonexistent", verbose = FALSE))
})

test_that("pivot_VDJ respects sample_col", {
  vdj <- generate_test_VDJ(n_cells = 20, n_clones = 5,
                            samples = c("S1", "S2"), seed = 104)
  wide <- pivot_VDJ(vdj, sample_col = "origin", verbose = FALSE)

  expect_true(all(wide$Sample %in% c("S1", "S2")))
})

test_that("pivot_VDJ TRA/TRB columns are V_CDR3_J or NA", {
  vdj <- generate_test_VDJ(n_cells = 30, n_clones = 5, seed = 105)
  wide <- pivot_VDJ(vdj, verbose = FALSE)

  # Non-NA values should contain underscores (V_CDR3_J format)
  for (col in c("TRA_1", "TRB_1")) {
    non_na <- wide[[col]][!is.na(wide[[col]])]
    if (length(non_na) > 0) {
      expect_true(all(grepl("_", non_na)),
                  info = paste(col, "should be V_CDR3_J format"))
    }
  }
})
