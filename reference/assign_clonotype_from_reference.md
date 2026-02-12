# Assign clonotypes by mapping to a reference table

Maps TCR data to an existing reference clonotype table using similarity
matrices at a specified stringency level. Useful for mapping new data
against previously identified clonotypes.

## Usage

``` r
assign_clonotype_from_reference(
  TCR_data,
  ref_table,
  similarity_matrix_type = .mapping.levels,
  clone_definition_df = .clone.definition.df
)
```

## Arguments

- TCR_data:

  A data frame with columns Sample, TRA_1, TRA_2, TRB_1, TRB_2.

- ref_table:

  A reference table with columns Sample, TRA_1, TRA_2, TRB_1, TRB_2, and
  CloneID (as output by
  [`identify_clonotypes`](https://doRemy.github.io/VDJigsaw/reference/identify_clonotypes.md)).

- similarity_matrix_type:

  Which stringency level to use for matching. One of the values in
  `.mapping.levels`.

- clone_definition_df:

  Clone definition data frame.

## Value

A character vector of mapped clone IDs, with `NA` for unmapped cells and
pipe-separated IDs for cells matching multiple clonotypes.
