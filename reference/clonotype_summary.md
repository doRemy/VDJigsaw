# Summarize clonotype changes across stringency levels

Tracks how clonotype assignments change as the stringency level is
relaxed, identifying which clonotypes are merged at each successive
level.

## Usage

``` r
clonotype_summary(TCR_data, mapping_levels = .mapping.levels)
```

## Arguments

- TCR_data:

  A data frame with clone ID columns for each stringency level (as
  output by
  [`assign_clonotype`](https://doRemy.github.io/VDJigsaw/reference/assign_clonotype.md)).

- mapping_levels:

  Character vector of stringency level names in order from most to least
  strict.

## Value

A data frame summarizing the number of clonotypes and whether merging
occurred at each stringency level transition.
