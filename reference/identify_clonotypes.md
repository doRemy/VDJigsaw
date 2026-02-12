# Identify clonotypes at multiple stringency levels

Groups cells into clonotypes using a hierarchy of stringency levels,
from strict (all four chains must match) to permissive (a single allele
suffices). Uses graph-based clustering via igraph to find connected
components of matching cells. Supports parallel processing across
samples.

## Usage

``` r
identify_clonotypes(
  TCR_data,
  clone_definition_df = .clone.definition.df,
  num_cores = 1,
  sample_col = NULL,
  clone_loose = "single_chain_single_allele",
  verbose = TRUE
)
```

## Arguments

- TCR_data:

  A data frame with columns TRA_1, TRA_2, TRB_1, TRB_2 containing
  validated TCR chain identifiers.

- clone_definition_df:

  Clone definition data frame specifying the stringency hierarchy.
  Defaults to the built-in definition.

- num_cores:

  Number of cores for parallel processing across samples. Set to `-1` to
  use all available cores. The actual number of cores used is capped at
  the number of available cores and the number of unique samples.

- sample_col:

  Name of the column to use as sample identifier. If `NULL`, all cells
  are treated as one sample.

- clone_loose:

  Which stringency level to use as the "loose" default. One of:
  "dual_chain_dual_allele", "dual_chain_one_partial",
  "dual_chain_both_partial", "single_chain_dual_allele",
  "single_chain_single_allele".

- verbose:

  Logical. If `TRUE`, print progress messages.

## Value

A list with two elements:

- TCR_data:

  Data frame with clone ID columns for each stringency level
  (CloneID.dual_chain_dual_allele, etc.) plus CloneID.loose.

- ref_tables:

  Named list of reference tables for each stringency level, each
  containing CloneID, Sample, and unique TCR chain info.
