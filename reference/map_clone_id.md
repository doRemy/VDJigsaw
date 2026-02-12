# Map clone IDs from VDJigsaw results to wide-format data

Extracts the clone ID columns from VDJigsaw results and aligns them to
the rows of the wide-format VDJ data frame.

## Usage

``` r
map_clone_id(
  VDJigsaw_res,
  VDJ_contigs_wide,
  cols_to_map = c("CloneID.dual_chain_dual_allele", "CloneID.loose",
    "CloneID.dual_chain_one_partial", "CloneID.dual_chain_both_partial",
    "CloneID.single_chain_dual_allele", "CloneID.single_chain_single_allele")
)
```

## Arguments

- VDJigsaw_res:

  Output from
  [`identify_clonotypes`](https://doRemy.github.io/VDJigsaw/reference/identify_clonotypes.md).

- VDJ_contigs_wide:

  Wide-format data frame (output of
  [`pivot_VDJ`](https://doRemy.github.io/VDJigsaw/reference/pivot_VDJ.md)).

- cols_to_map:

  Character vector of clone ID column names to extract.

## Value

A data frame with the requested clone ID columns, aligned to the rows of
`VDJ_contigs_wide`.
