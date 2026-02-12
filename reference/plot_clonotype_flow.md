# Plot clonotype flow across stringency levels

Creates an alluvial diagram showing how clonotypes merge as stringency
is relaxed. Each vertical axis represents a stringency level, strata
represent clonotypes, and flows show cells moving between clone
definitions. Top-N filtering keeps the diagram readable.

## Usage

``` r
plot_clonotype_flow(
  TCR_data,
  mapping_levels = .mapping.levels,
  top_n = 20,
  sample_filter = NULL,
  per_sample = FALSE,
  combined = TRUE,
  show_unassigned = FALSE,
  title = "Clonotype Flow Across Stringency Levels"
)
```

## Arguments

- TCR_data:

  A data frame with clone ID columns for each stringency level (as
  output by
  [`assign_clonotype`](https://doRemy.github.io/VDJigsaw/reference/assign_clonotype.md)).

- mapping_levels:

  Character vector of stringency level names in order from most to least
  strict.

- top_n:

  Maximum number of clones to show at the loosest level. All upstream
  clones feeding into these are kept; the rest are grouped as "Other".

- sample_filter:

  Optional sample name to filter to. If `NULL`, all samples are
  included. Ignored when `per_sample = TRUE`.

- per_sample:

  Logical. If `TRUE`, create a separate plot for each unique sample.

- combined:

  Logical. If `TRUE` and `per_sample = TRUE`, combine all per-sample
  plots into a single figure using `patchwork`. If `FALSE`, return a
  named list of plots.

- show_unassigned:

  Logical. If `TRUE`, show "Unassigned" strata for cells with `NA` clone
  IDs at stricter levels.

- title:

  Plot title.

## Value

A `ggplot` object, or a named list of `ggplot` objects when
`per_sample = TRUE` and `combined = FALSE`.

## Details

Requires the `ggalluvial` package.
