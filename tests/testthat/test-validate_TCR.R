test_that("validate_TCR adds .invalid columns", {
  vdj <- generate_test_VDJ(n_cells = 20, n_clones = 5, seed = 200)
  wide <- pivot_VDJ(vdj, verbose = FALSE)
  validated <- validate_TCR(wide, verbose = FALSE)

  expect_true("TRA_1.invalid" %in% colnames(validated))
  expect_true("TRA_2.invalid" %in% colnames(validated))
  expect_true("TRB_1.invalid" %in% colnames(validated))
  expect_true("TRB_2.invalid" %in% colnames(validated))
  expect_true(all(is.logical(validated$TRA_1.invalid)))
})

test_that("validate_TCR removes duplicate alleles", {
  # Use generated data and then manually create duplicates
  vdj <- generate_test_VDJ(n_cells = 20, n_clones = 5, seed = 202)
  wide <- pivot_VDJ(vdj, verbose = FALSE)

  # Force a duplicate: set TRA_2 = TRA_1 for first non-NA row
  idx <- which(!is.na(wide$TRA_1))[1]
  wide$TRA_2[idx] <- wide$TRA_1[idx]

  result <- validate_TCR(wide, verbose = FALSE)
  # TRA_2 should be NA after validation (duplicate removed)
  expect_true(is.na(result$TRA_2[idx]))
})

test_that("validate_TCR swaps alleles when allele_1 is NA but allele_2 is not", {
  vdj <- generate_test_VDJ(n_cells = 20, n_clones = 5, seed = 203)
  wide <- pivot_VDJ(vdj, verbose = FALSE)

  # Force a swap scenario: set TRA_1 = NA and TRA_2 = some value
  idx <- which(!is.na(wide$TRA_1))[1]
  original_value <- wide$TRA_1[idx]
  wide$TRA_2[idx] <- original_value
  wide$TRA_1[idx] <- NA

  result <- validate_TCR(wide, verbose = FALSE)
  # TRA_1 should now have the value, TRA_2 should be NA
  expect_equal(result$TRA_1[idx], original_value)
  expect_true(is.na(result$TRA_2[idx]))
})

test_that("validate_TCR errors on missing columns", {
  df <- data.frame(TRA_1 = "x", TRA_2 = "y")
  expect_error(validate_TCR(df, verbose = FALSE), "Missing required columns")
})

test_that("validate_TCR handles remove_invalid_VDJ option", {
  vdj <- generate_test_VDJ(n_cells = 50, n_clones = 10, invalid_rate = 0.3, seed = 201)
  wide <- pivot_VDJ(vdj, verbose = FALSE)

  result_keep <- validate_TCR(wide, remove_invalid_VDJ = FALSE, verbose = FALSE)
  result_remove <- validate_TCR(wide, remove_invalid_VDJ = TRUE, verbose = FALSE)

  # With removal, should have more NAs (or equal) than without
  na_keep <- sum(is.na(result_keep$TRA_1)) + sum(is.na(result_keep$TRB_1))
  na_remove <- sum(is.na(result_remove$TRA_1)) + sum(is.na(result_remove$TRB_1))
  expect_true(na_remove >= na_keep)
})
