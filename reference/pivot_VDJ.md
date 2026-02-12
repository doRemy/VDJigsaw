# Pivot VDJ contig data to wide format

Converts VDJ contig annotation data from long format (one row per
contig) to wide format (one row per cell barcode), creating composite
TCR chain identifiers (TRA_1, TRA_2, TRB_1, TRB_2) from V gene, CDR3,
and J gene.

## Usage

``` r
pivot_VDJ(
  VDJ_data,
  sample_col = NULL,
  is_cell = FALSE,
  high_confidence = FALSE,
  productive = FALSE,
  full_length = FALSE,
  cdr3_col = "cdr3_nt",
  verbose = TRUE
)
```

## Arguments

- VDJ_data:

  A data frame of VDJ contig annotations in 10x Genomics format. Must
  contain columns: barcode, chain, v_gene, cdr3, j_gene, umis, and other
  standard 10x VDJ columns.

- sample_col:

  Name of the column in `VDJ_data` to use as sample identifier. If
  `NULL`, all cells are assigned to a single sample called "Sample".

- is_cell:

  Logical. If `TRUE`, filter to rows where `is_cell == "true"`.

- high_confidence:

  Logical. If `TRUE`, filter to rows where `high_confidence == "true"`.

- productive:

  Logical. If `TRUE`, filter to rows where `productive == "true"`.

- full_length:

  Logical. If `TRUE`, filter to rows where `full_length == "true"`.

- cdr3_col:

  Which CDR3 column to use for building the composite TCR chain
  identifier. Either `"cdr3_nt"` (nucleotide sequence, default) or
  `"cdr3"` (amino acid sequence).

- verbose:

  Logical. If `TRUE`, print progress messages.

## Value

A data frame in wide format with one row per cell barcode. Contains
columns for sample, barcode, TCR chain identifiers (TRA_1, TRA_2, TRB_1,
TRB_2), contig metadata, and individual chain component columns.
