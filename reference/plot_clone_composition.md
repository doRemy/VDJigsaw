# Plot clonotype chain composition

Creates a tile chart showing the V gene usage for each chain position
(TRA_1, TRA_2, TRB_1, TRB_2) across clonotypes at a given stringency
level. Each tile is colored by V gene identity, making it easy to see
which clones share V gene segments and which chains are present or
absent.

## Usage

``` r
plot_clone_composition(
  ref_table,
  top_n = 20,
  sample_filter = NULL,
  per_sample = FALSE,
  combined = TRUE,
  title = "Clonotype Chain Composition"
)
```

## Arguments

- ref_table:

  A reference table from
  [`identify_clonotypes`](https://doRemy.github.io/VDJigsaw/reference/identify_clonotypes.md)
  or
  [`assign_clonotype`](https://doRemy.github.io/VDJigsaw/reference/assign_clonotype.md)
  output (e.g., `result$ref_tables$dual_chain_both_partial`). Must
  contain columns: TRA_1, TRB_1.

- top_n:

  Maximum number of clones to display. Clones are ranked by chain
  richness (number of non-NA chains). Remaining clones are excluded.

- sample_filter:

  Optional sample name to filter to. If `NULL`, all samples are
  included. Ignored when `per_sample = TRUE`.

- per_sample:

  Logical. If `TRUE`, create a separate plot for each unique sample in
  the reference table.

- combined:

  Logical. If `TRUE` and `per_sample = TRUE`, combine all per-sample
  plots into a single figure using `patchwork`. If `FALSE`, return a
  named list of plots.

- title:

  Plot title.

## Value

A `ggplot` object, or a named list of `ggplot` objects when
`per_sample = TRUE` and `combined = FALSE`.
