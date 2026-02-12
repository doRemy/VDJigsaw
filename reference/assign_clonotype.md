# Assign clonotypes from raw VDJ contig data

Complete VDJigsaw pipeline that takes raw VDJ contig annotations, pivots
them to wide format, validates TCR chains, identifies clonotypes at
multiple stringency levels, and returns the annotated data with clone
assignments.

## Usage

``` r
assign_clonotype(
  VDJ_data,
  sample_col = NULL,
  clone_definition_df = .clone.definition.df,
  is_cell = FALSE,
  high_confidence = FALSE,
  productive = FALSE,
  full_length = FALSE,
  remove_invalid_VDJ = FALSE,
  cdr3_col = "cdr3_nt",
  num_cores = 1,
  clone_loose = "single_chain_single_allele",
  verbose = TRUE
)
```

## Arguments

- VDJ_data:

  A data frame of VDJ contig annotations in 10x Genomics format.

- sample_col:

  Name of the column in `VDJ_data` to use as sample identifier. If
  `NULL`, all cells are assigned to "Sample".

- clone_definition_df:

  Clone definition data frame specifying the stringency hierarchy.
  Defaults to the built-in definition.

- is_cell:

  Logical. If `TRUE`, filter to rows where `is_cell == "true"`.

- high_confidence:

  Logical. If `TRUE`, filter to rows where `high_confidence == "true"`.

- productive:

  Logical. If `TRUE`, filter to rows where `productive == "true"`.

- full_length:

  Logical. If `TRUE`, filter to rows where `full_length == "true"`.

- remove_invalid_VDJ:

  Logical. If `TRUE`, set invalid VDJ sequences to `NA`.

- cdr3_col:

  Which CDR3 column to use for building the composite TCR chain
  identifier. Either `"cdr3_nt"` (nucleotide, default) or `"cdr3"`
  (amino acid).

- num_cores:

  Number of cores for parallel processing. Set to `-1` to use all
  available cores. Capped at the number of available cores and the
  number of unique samples.

- clone_loose:

  Which stringency level to use as the default "loose" clone definition.

- verbose:

  Logical. If `TRUE`, print progress messages.

## Value

A list with two elements:

- TCR_data:

  Wide-format data frame with TCR chains, clone IDs at each stringency
  level, and a default CloneID column.

- ref_tables:

  Named list of reference tables for each stringency level.

## Examples

``` r
if (FALSE) { # \dontrun{
VDJ_contigs <- read.csv("filtered_contig_annotations.csv")
result <- assign_clonotype(VDJ_contigs, sample_col = "origin")
TCR_data <- result$TCR_data
ref_tables <- result$ref_tables
} # }
```
