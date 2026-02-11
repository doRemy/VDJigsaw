test_that("assign_clonotype_from_reference maps cells to existing clones", {
  vdj <- generate_test_VDJ(n_cells = 30, n_clones = 5, seed = 600)

  # Build the query data with pivot + validate (has Sample, TRA_1, etc.)
  wide <- pivot_VDJ(vdj, sample_col = "origin", verbose = FALSE)
  validated <- validate_TCR(wide, verbose = FALSE)

  # Run clonotype identification to get reference tables
  clono_res <- identify_clonotypes(validated, sample_col = "Sample",
                                    num_cores = 1, verbose = FALSE)
  ref <- clono_res$ref_tables$dual_chain_dual_allele

  # Map validated data (which has Sample, TRA_1, TRA_2, TRB_1, TRB_2)
  mapped <- assign_clonotype_from_reference(
    TCR_data = validated,
    ref_table = ref,
    similarity_matrix_type = "dual_chain_dual_allele"
  )

  expect_type(mapped, "character")
  expect_equal(length(mapped), nrow(validated))
})

test_that("assign_clonotype_from_reference errors on missing columns", {
  df <- data.frame(TRA_1 = "x")
  ref <- data.frame(TRA_1 = "x")
  expect_error(assign_clonotype_from_reference(df, ref), "Missing required")
})

test_that("annotate_by_clone returns correct columns", {
  vdj <- generate_test_VDJ(n_cells = 30, n_clones = 5, seed = 601)
  result <- assign_clonotype(vdj, num_cores = 1, verbose = FALSE)

  ref <- result$ref_tables$single_chain_single_allele
  clone_ids <- result$TCR_data$CloneID.single_chain_single_allele
  clone_ids <- clone_ids[!is.na(clone_ids)]

  if (length(clone_ids) > 0) {
    annotated <- annotate_by_clone(clone_ids, ref)

    expect_s3_class(annotated, "data.frame")
    expect_equal(ncol(annotated), ncol(ref))
    expect_equal(nrow(annotated), length(clone_ids))
  }
})
