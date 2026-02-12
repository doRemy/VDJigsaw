# Plot clonotype assignment summary across stringency levels

Creates a two-panel bar chart showing how clonotype counts and cell
assignments change across stringency levels. The top panel shows the
number of unique clonotypes at each level (decreasing as clones merge).
The bottom panel shows the number of cells assigned vs unassigned at
each level. Counts are displayed on top of each bar.

## Usage

``` r
plot_stringency_summary(
  TCR_data,
  mapping_levels = .mapping.levels,
  sample_filter = NULL,
  per_sample = FALSE,
  combined = TRUE,
  title = "Clonotype Assignment Summary"
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

- sample_filter:

  Optional sample name to filter to. If `NULL`, all samples are
  included. Ignored when `per_sample = TRUE`.

- per_sample:

  Logical. If `TRUE`, create a separate plot for each unique sample.

- combined:

  Logical. If `TRUE` and `per_sample = TRUE`, combine all per-sample
  plots into a single figure using `patchwork`. If `FALSE`, return a
  named list of plots.

- title:

  Plot title.

## Value

A `ggplot` object, or a named list of `ggplot` objects when
`per_sample = TRUE` and `combined = FALSE`.
