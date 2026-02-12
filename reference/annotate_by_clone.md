# Annotate clone IDs with reference table metadata

Given a vector of clone IDs (possibly pipe-separated for multi-mapped
cells), retrieves corresponding annotations from a reference table.
Handles numeric, logical, factor, and character columns appropriately.

## Usage

``` r
annotate_by_clone(clone_ids, ref_table)
```

## Arguments

- clone_ids:

  Character vector of clone IDs. Pipe-separated IDs indicate cells
  mapping to multiple clonotypes.

- ref_table:

  A reference table with clone IDs as row names and annotation columns.

## Value

A data frame with the same columns as `ref_table`, one row per input
clone ID.
