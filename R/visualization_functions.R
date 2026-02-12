# Visualization functions for VDJigsaw
# @Author Rémy Pétremand

# ------------------------------------------------------------------------------
# Standard function to check presence of essential libraries: ggplot2, ggalluvial and patchwork
# ------------------------------------------------------------------------------

.check_ggplot2_available <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function. ",
         "Install it with: install.packages('ggplot2')",
         call. = FALSE)
  }
}
.check_ggalluvial_available <- function() {
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    stop("Package 'ggalluvial' is required for this function. ",
         "Install it with: install.packages('ggalluvial')",
         call. = FALSE)
  }
}
.check_patchwork_available <- function() {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required when combined = TRUE. ",
         "Install it with: install.packages('patchwork')",
         call. = FALSE)
  }
}

# ------------------------------------------------------------------------------
# Summary table
# ------------------------------------------------------------------------------

#' Summarize clonotype changes across stringency levels
#'
#' Tracks how clonotype assignments change as the stringency level is relaxed,
#' identifying which clonotypes are merged at each successive level.
#'
#' @param TCR_data A data frame with clone ID columns for each stringency
#'   level (as output by \code{\link{assign_clonotype}}).
#' @param mapping_levels Character vector of stringency level names in order
#'   from most to least strict.
#'
#' @return A data frame summarizing the number of clonotypes and whether
#'   merging occurred at each stringency level transition.
#'
#' @export
clonotype_summary <- function(TCR_data, mapping_levels = .mapping.levels){
  
  # Part 1: Get the vector of columns and check that everything is in order
  clone.columns <- paste0("CloneID.", mapping_levels)
  track.change.columns <- paste0("Match.", mapping_levels)
  
  if (!all(clone.columns %in% colnames(TCR_data))){
    missing.cols <- clone.columns[!clone.columns %in% colnames(TCR_data)]
    stop("Missing required clone ID columns in TCR_data: ", paste(missing.cols, collapse = ", "),
         ". Available columns: ", paste(colnames(TCR_data), collapse = ", "))
  }
  
  # Part 2: Set up the TCR_data data.frame with the track.change.columns
  TCR_data <- TCR_data %>% 
    group_by_at(clone.columns) %>%
    summarise(n = n(), .groups = "keep") %>% 
    as.data.frame()
  TCR_data[,track.change.columns] <- matrix(FALSE, nrow = nrow(TCR_data), ncol = length(track.change.columns))
  TCR_data[,track.change.columns[1]] <- TRUE
  
  # Part 3: Track change
  for (i in 2:length(clone.columns)){
    count.df <- TCR_data %>%
      group_by_at(clone.columns[i]) %>%
      summarise(n_dup = length(unique(.data[[clone.columns[i-1]]])), .groups = "keep") %>%
      filter(n_dup > 1) %>%
      as.data.frame()
    mask.dup <- TCR_data[[clone.columns[i]]] %in% count.df[[clone.columns[i]]]
    TCR_data[mask.dup, track.change.columns[i]] <- TRUE
  }
  
  # Part 4: Get the summary and return it
  summary.df <- TCR_data %>% 
    group_by_at(track.change.columns) %>% 
    summarise(n_clones = n(), .groups = "keep") %>% 
    as.data.frame()
  
  return(summary.df)
}

# ------------------------------------------------------------------------------
# Stringency summary plot
# ------------------------------------------------------------------------------

#' Plot clonotype assignment summary across stringency levels
#'
#' Creates a two-panel bar chart showing how clonotype counts and cell
#' assignments change across stringency levels. The top panel shows the
#' number of unique clonotypes at each level (decreasing as clones merge).
#' The bottom panel shows the number of cells assigned vs unassigned at
#' each level. Counts are displayed on top of each bar.
#'
#' @param TCR_data A data frame with clone ID columns for each stringency
#'   level (as output by \code{\link{assign_clonotype}}).
#' @param mapping_levels Character vector of stringency level names in order
#'   from most to least strict.
#' @param sample_filter Optional sample name to filter to. If \code{NULL},
#'   all samples are included. Ignored when \code{per_sample = TRUE}.
#' @param per_sample Logical. If \code{TRUE}, create a separate plot for
#'   each unique sample.
#' @param combined Logical. If \code{TRUE} and \code{per_sample = TRUE},
#'   combine all per-sample plots into a single figure using
#'   \code{patchwork}. If \code{FALSE}, return a named list of plots.
#' @param title Plot title.
#'
#' @return A \code{ggplot} object, or a named list of \code{ggplot} objects
#'   when \code{per_sample = TRUE} and \code{combined = FALSE}.
#'
#' @export
plot_stringency_summary <- function(
    TCR_data,
    mapping_levels = .mapping.levels,
    sample_filter = NULL,
    per_sample = FALSE,
    combined = TRUE,
    title = "Clonotype Assignment Summary") {
  
  .check_ggplot2_available()
  
  # Validate columns
  clone_columns <- paste0("CloneID.", mapping_levels)
  if (!all(clone_columns %in% colnames(TCR_data))) {
    missing_cols <- clone_columns[!clone_columns %in% colnames(TCR_data)]
    stop("Missing required clone ID columns in TCR_data: ", paste(missing_cols, collapse = ", "))
  }
  
  # Per-sample mode
  if (per_sample) {
    if (!"Sample" %in% colnames(TCR_data)) {
      stop("Sample column not found in TCR_data. Cannot split by sample.")
    }
    samples <- sort(unique(TCR_data$Sample))
    plot_list <- lapply(setNames(samples, samples), function(s) {
      plot_stringency_summary(
        TCR_data = TCR_data, 
        mapping_levels = mapping_levels,
        sample_filter = s, 
        per_sample = FALSE, 
        combined = FALSE,
        title = paste0(title, " - ", s))
    })
    if (combined) {
      .check_patchwork_available()
      return(patchwork::wrap_plots(plot_list, ncol = min(length(plot_list), 2)))
    }
    return(plot_list)
  }
  
  # Filter by sample if requested
  if (!is.null(sample_filter)) {
    if (!"Sample" %in% colnames(TCR_data)) {
      stop("Sample column not found in TCR_data")
    }
    TCR_data <- TCR_data[TCR_data$Sample == sample_filter, ]
    if (nrow(TCR_data) == 0) {
      stop("No data remaining after filtering to sample '", sample_filter, "'")
    }
  }
  
  n_total <- nrow(TCR_data)
  
  # Compute summary stats per level
  summary_df <- TCR_data %>%
    select(all_of(clone_columns)) %>%
    pivot_longer(cols = all_of(clone_columns), names_to = "clone_col", values_to = "clone_id") %>%
    group_by(clone_col) %>%
    summarise(
      n_clones     = n_distinct(clone_id, na.rm = TRUE),
      n_assigned   = sum(!is.na(clone_id)),
      n_unassigned = sum(is.na(clone_id)),
      .groups = "drop"
    ) %>%
    mutate(
      level_label  = factor(
        gsub("_", " ", gsub("^CloneID\\.", "", clone_col)),
        levels = gsub("_", " ", mapping_levels)),
      pct_assigned = round(100 * n_assigned / n_total, 1)
    )
  
  # Build long-format data for faceted plot
  plot_df <- bind_rows(
    summary_df %>%
      transmute(
        level_label = level_label,
        metric = factor("Unique Clones", levels = c("Unique Clones", "Cells Assigned")),
        value = n_clones,
        label = as.character(n_clones)
      ),
    summary_df %>%
      transmute(
        level_label = level_label,
        metric = factor("Cells Assigned", levels = c("Unique Clones", "Cells Assigned")),
        value  = n_assigned,
        label  = paste0(n_assigned, " (", pct_assigned, "%)")
      )
  )
  
  # Create plot
  p <- ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(x = level_label, y = value)) +
    ggplot2::geom_col(fill = "#008c93") +
    ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.5, size = 3) +
    ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 1) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(title = title, x = "Stringency level", y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1, size = 9),
      strip.text = ggplot2::element_text(face = "bold", size = 11)
    )

  return(p)
}

# ------------------------------------------------------------------------------
# Clone composition plot
# ------------------------------------------------------------------------------

#' Plot clonotype chain composition
#'
#' Creates a tile chart showing the V gene usage for each chain position
#' (TRA_1, TRA_2, TRB_1, TRB_2) across clonotypes at a given stringency
#' level. Each tile is colored by V gene identity, making it easy to see
#' which clones share V gene segments and which chains are present or absent.
#'
#' @param ref_table A reference table from
#'   \code{\link{identify_clonotypes}} or \code{\link{assign_clonotype}}
#'   output (e.g., \code{result$ref_tables$dual_chain_both_partial}).
#'   Must contain columns: TRA_1, TRB_1.
#' @param top_n Maximum number of clones to display. Clones are ranked by
#'   chain richness (number of non-NA chains). Remaining clones are excluded.
#' @param sample_filter Optional sample name to filter to. If \code{NULL},
#'   all samples are included. Ignored when \code{per_sample = TRUE}.
#' @param per_sample Logical. If \code{TRUE}, create a separate plot for
#'   each unique sample in the reference table.
#' @param combined Logical. If \code{TRUE} and \code{per_sample = TRUE},
#'   combine all per-sample plots into a single figure using
#'   \code{patchwork}. If \code{FALSE}, return a named list of plots.
#' @param title Plot title.
#'
#' @return A \code{ggplot} object, or a named list of \code{ggplot} objects
#'   when \code{per_sample = TRUE} and \code{combined = FALSE}.
#'
#' @export
plot_clone_composition <- function(
    ref_table,
    top_n = 20,
    sample_filter = NULL,
    per_sample = FALSE,
    combined = TRUE,
    title = "Clonotype Chain Composition") {
  
  .check_ggplot2_available()
  
  # Validate columns
  required_cols <- c("TRA_1", "TRB_1")
  if (!all(required_cols %in% colnames(ref_table))) {
    missing_cols <- required_cols[!required_cols %in% colnames(ref_table)]
    stop("Missing required columns in ref_table: ", paste(missing_cols, collapse = ", "))
  }
  
  # Ensure CloneID column exists (it's the rownames in ref_tables)
  if (!"CloneID" %in% colnames(ref_table)) {
    ref_table$CloneID <- rownames(ref_table)
  }
  
  # Per-sample mode
  if (per_sample) {
    if (!"Sample" %in% colnames(ref_table)) {
      stop("Sample column not found in ref_table. Cannot split by sample.")
    }
    samples <- sort(unique(ref_table$Sample))
    plot_list <- lapply(setNames(samples, samples), function(s) {
      plot_clone_composition(
        ref_table, 
        top_n = top_n,
        sample_filter = s, 
        per_sample = FALSE, 
        combined = FALSE,
        title = paste0(title, " - ", s)
        )
    })
    if (combined) {
      .check_patchwork_available()
      return(patchwork::wrap_plots(plot_list, ncol = min(length(plot_list), 2)))
    }
    return(plot_list)
  }
  
  # Filter by sample if requested
  if (!is.null(sample_filter)) {
    if (!"Sample" %in% colnames(ref_table)) {
      stop("Sample column not found in ref_table")
    }
    ref_table <- ref_table[ref_table$Sample == sample_filter, ]
    if (nrow(ref_table) == 0) {
      stop("No data remaining after filtering to sample '", sample_filter, "'")
    }
  }
  
  # Compute chain richness, rank, and keep top_n
  chain_cols <- intersect(c("TRA_1", "TRA_2", "TRB_1", "TRB_2"), colnames(ref_table))
  ref_table <- ref_table %>%
    mutate(n_chains = rowSums(!is.na(across(all_of(chain_cols))))) %>%
    arrange(desc(n_chains)) %>%
    slice(seq_len(min(top_n, n())))
  
  # Reshape to long format and extract V gene
  ref_long <- ref_table %>%
    select("CloneID", all_of(chain_cols)) %>%
    pivot_longer(cols = all_of(chain_cols), names_to = "Chain", values_to = "Value") %>%
    mutate(
      V_gene  = ifelse(is.na(Value), NA_character_, sub("_.*", "", Value)),
      Present = !is.na(Value),
      CloneID = factor(CloneID, levels = rev(ref_table$CloneID)),
      Chain   = factor(Chain, levels = chain_cols)
    )
  
  # Create Plot
  p <- ggplot2::ggplot(ref_long, ggplot2::aes(x = Chain, y = CloneID)) +
    ggplot2::geom_tile(ggplot2::aes(fill = V_gene), color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(Present, V_gene, "")), size = 2, na.rm = TRUE) +
    ggplot2::scale_fill_discrete(na.value = "grey95") +
    ggplot2::labs(title = title, x = "Chain position", y = "Clone ID", fill = "V gene") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 7)
    )

  return(p)
}

# ------------------------------------------------------------------------------
# Clonotype flow (alluvial) plot
# ------------------------------------------------------------------------------

#' Plot clonotype flow across stringency levels
#'
#' Creates an alluvial diagram showing how clonotypes merge as stringency
#' is relaxed. Each vertical axis represents a stringency level, strata
#' represent clonotypes, and flows show cells moving between clone
#' definitions. Top-N filtering keeps the diagram readable.
#'
#' Requires the \code{ggalluvial} package.
#'
#' @param TCR_data A data frame with clone ID columns for each stringency
#'   level (as output by \code{\link{assign_clonotype}}).
#' @param mapping_levels Character vector of stringency level names in order
#'   from most to least strict.
#' @param top_n Maximum number of clones to show at the loosest level. All
#'   upstream clones feeding into these are kept; the rest are grouped as
#'   "Other".
#' @param sample_filter Optional sample name to filter to. If \code{NULL},
#'   all samples are included. Ignored when \code{per_sample = TRUE}.
#' @param per_sample Logical. If \code{TRUE}, create a separate plot for
#'   each unique sample.
#' @param combined Logical. If \code{TRUE} and \code{per_sample = TRUE},
#'   combine all per-sample plots into a single figure using
#'   \code{patchwork}. If \code{FALSE}, return a named list of plots.
#' @param show_unassigned Logical. If \code{TRUE}, show "Unassigned" strata
#'   for cells with \code{NA} clone IDs at stricter levels.
#' @param title Plot title.
#'
#' @return A \code{ggplot} object, or a named list of \code{ggplot} objects
#'   when \code{per_sample = TRUE} and \code{combined = FALSE}.
#'
#' @export
plot_clonotype_flow <- function(
    TCR_data,
    mapping_levels = .mapping.levels,
    top_n = 20,
    sample_filter = NULL,
    per_sample = FALSE,
    combined = TRUE,
    show_unassigned = FALSE,
    title = "Clonotype Flow Across Stringency Levels") {
  
  .check_ggplot2_available()
  .check_ggalluvial_available()
  
  # Validate columns
  clone_columns <- paste0("CloneID.", mapping_levels)
  if (!all(clone_columns %in% colnames(TCR_data))) {
    missing_cols <- clone_columns[!clone_columns %in% colnames(TCR_data)]
    stop("Missing required clone ID columns in TCR_data: ", paste(missing_cols, collapse = ", "))
  }
  
  # Per-sample mode
  if (per_sample) {
    if (!"Sample" %in% colnames(TCR_data)) {
      stop("Sample column not found in TCR_data. Cannot split by sample.")
    }
    samples <- sort(unique(TCR_data$Sample))
    plot_list <- lapply(setNames(samples, samples), function(s) {
      plot_clonotype_flow(
        TCR_data = TCR_data, 
        mapping_levels = mapping_levels, 
        top_n = top_n,
        sample_filter = s, 
        per_sample = FALSE, 
        combined = FALSE,
        show_unassigned = show_unassigned,
        title = paste0(title, " - ", s))
    })
    if (combined) {
      .check_patchwork_available()
      return(patchwork::wrap_plots(plot_list, ncol = 1))
    }
    return(plot_list)
  }
  
  # Filter by sample
  if (!is.null(sample_filter)) {
    if (!"Sample" %in% colnames(TCR_data)) {
      stop("Sample column not found in TCR_data")
    }
    TCR_data <- TCR_data[TCR_data$Sample == sample_filter, ]
    if (nrow(TCR_data) == 0) {
      stop("No data remaining after filtering to sample '", sample_filter, "'")
    }
  }
  
  # Work with clone columns only and handle NAs
  work_data <- TCR_data[, clone_columns, drop = F]
  if (!show_unassigned) {
    work_data <- work_data %>%
      filter(rowSums(!is.na(across(all_of(clone_columns)))) > 0)
  }
  work_data <- work_data %>%
    mutate(across(all_of(clone_columns), ~ if_else(is.na(.), "Unassigned", .)))

  if (nrow(work_data) == 0) {
    stop("No cells with clonotype assignments found")
  }

  # Top-N filtering at loosest level
  # Get the loosest clone definition column as the last element of clone_columns
  loosest_col <- clone_columns[length(clone_columns)]
  keep_clones <- work_data %>%
    count(.data[[loosest_col]], name = "n") %>%
    filter(.data[[loosest_col]] != "Unassigned") %>%
    arrange(desc(.data$n)) %>%
    slice(seq_len(min(top_n, n()))) %>%
    pull(.data[[loosest_col]])

  work_data <- work_data %>%
    mutate(across(all_of(clone_columns),
          ~ if_else(
              condition = .data[[loosest_col]] %in% keep_clones | .data[[loosest_col]] == "Unassigned",
              true = .,
              false = "Other"
            )
          )
        )

  # Aggregate: count cells per unique combination of clone assignments
  agg <- work_data %>%
    group_by(across(all_of(clone_columns))) %>%
    summarise(Freq = n(), .groups = "drop") %>%
    mutate(
      fill_clone  = .data[[loosest_col]],
      alluvium_id = seq_len(n())
    )

  # Pivot to long format for ggalluvial
  level_labels <- gsub("_", " ", gsub("^CloneID\\.", "", clone_columns))

  agg_long <- agg %>%
    pivot_longer(
      cols = all_of(clone_columns),
      names_to = "level",
      values_to = "stratum"
    ) %>%
    mutate(level = factor(level, levels = clone_columns, labels = level_labels))

  # Get the plot
  p <- ggplot2::ggplot(
    data = agg_long,
    mapping = ggplot2::aes(
      x = level, 
      y = Freq, 
      stratum = stratum, 
      alluvium = alluvium_id, 
      fill = fill_clone)) +
    ggalluvial::geom_flow(alpha = 0.6) +
    ggalluvial::geom_stratum(width = 1/3, color = "grey20", linewidth = 0.2) +
    ggplot2::labs(title = title, x = "Stringency level", y = "Number of cells", fill = "Clone (loosest)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1, size = 9),
      legend.position = "right"
    )

  return(p)
}