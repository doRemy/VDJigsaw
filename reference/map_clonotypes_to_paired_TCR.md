# Map clonotypes from allele format to paired TCR format

Maps TCR data from individual allele format (TRA_1/TRA_2/TRB_1/TRB_2) to
paired alpha/beta chain format (TRA/TRB). Useful when the reference
table uses only TRA and TRB (no second allele) and you want to map
allele-level data to this paired format.

## Usage

``` r
map_clonotypes_to_paired_TCR(
  data_to_map,
  ref_table,
  col_to_map,
  map_to_label = NULL
)
```

## Arguments

- data_to_map:

  A data frame with columns Sample, TRA_1, TRA_2, TRB_1, TRB_2.

- ref_table:

  A reference table with columns Sample, TRA, TRB, and the column
  specified by `col_to_map`.

- col_to_map:

  Name of the column in `ref_table` to map values from.

- map_to_label:

  Optional data frame for additional label mapping. Row names should
  correspond to mapped values. Factor and non-factor columns are handled
  differently.

## Value

A data frame with mapping results including strict and loose matches,
mapping type (strict/loose_A/loose_B/ND), and optionally label-mapped
columns.
