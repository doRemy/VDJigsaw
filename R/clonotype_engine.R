# Clonotype engine: internal helpers, identify_clonotypes, assign_clonotype, map_clone_id
# @Author Rémy Pétremand

# ------------------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------------------

#' Validate clone definition dataframe
#' @param clone_definition_df A data frame with clone definition columns.
#' @keywords internal
#' @noRd
.check_clone_definition_df <- function(clone_definition_df){
  necessary.cols <- colnames(.clone.definition.df)
  if (!all(necessary.cols %in% colnames(clone_definition_df))){
    missing.cols <- necessary.cols[!necessary.cols %in% colnames(clone_definition_df)]
    stop("Missing required columns in clone_definition_df: ", paste(missing.cols, collapse = ", "),
         ". Expected: ", paste(necessary.cols, collapse = ", "))
  }
  if (!all(unique(clone_definition_df$CloneDef) %in% .clone.definition.df$CloneDef)){
    warning("[VDJigsaw] Some clone definitions are not in the default set (might cause issues)")
  }
  test <- clone_definition_df %>% group_by(.data$CloneDef) %>% reframe(n = length(unique(.data$CloneDefRef)))
  if (any(test$n > 1)){
    stop("clone_definition_df error: Mismatched CloneDef and CloneDefRef. ",
         "Each CloneDef should have exactly one CloneDefRef")
  }
}

#' Extract unique TCR elements for a given chain
#' @param x1 First allele values.
#' @param x2 Second allele values.
#' @param TCR_chain Chain name prefix (e.g., "TRA" or "TRB").
#' @return A one-row tibble with chain count and sorted alleles.
#' @keywords internal
#' @noRd
.unique_TCR_element <- function(x1, x2, TCR_chain){
  # Remove NA values, sort alphabetically and ensure at least 2 elements
  x <- c(x1, x2)
  x <- c(sort(unique(x[!is.na(x)])), as.character(NA), as.character(NA))
  # Get results by selecting the two first chains
  res <- tibble(n_chain = (length(x) - 2), AB_1 = x[1], AB_2 = x[2])
  names(res) <- paste(TCR_chain, c("n_chain", "1", "2"), sep = "_")
  return(res)
}

#' Generate clone version columns for all stringency combinations
#' @param TCR_data Data frame with TRA_1, TRA_2, TRB_1, TRB_2, Sample columns.
#' @param clone_definition_df Clone definition data frame.
#' @return Data frame with clone version columns.
#' @keywords internal
#' @noRd
.get_clone_versions <- function(TCR_data, clone_definition_df = .clone.definition.df){

  # Sanity check:
  .check_clone_definition_df(clone_definition_df)

  # Count the number of TRA and TRB defined
  TCR_data$nTRA <- rowSums(!is.na(TCR_data[,c("TRA_1", "TRA_2")]))
  TCR_data$nTRB <- rowSums(!is.na(TCR_data[,c("TRB_1", "TRB_2")]))

  # Set up the NAcol:
  TCR_data$NAcol <- as.character(NA)

  # Get the clone versions:
  for (row_i in 1:nrow(clone_definition_df)){
    # Set the name of the columns of interest
    new.col <- clone_definition_df[row_i, "ColName"]
    TCR.cols <- as.character(clone_definition_df[row_i, c("TRA_1", "TRA_2", "TRB_1", "TRB_2")])
    clone.cols <- c("Sample", TCR.cols)

    # Get the clone versions
    TCR_data[[new.col]] <- apply(X = TCR_data[,clone.cols], MARGIN = 1, FUN = paste, sep = "_", collapse = "_")

    # Remove the aberrant clone versions 
    # 1. When clone version is <Sample>_NA_NA_NA_NA to NA
    TCR_data[[new.col]][rowSums(is.na(TCR_data[,c("TRA_1", "TRA_2", "TRB_1", "TRB_2")])) == 4] <- as.character(NA)
    # 2. When the second allele is defined and the first one not 
    mask.2.rm <- (is.na(TCR_data[,TCR.cols[1]]) & (!is.na(TCR_data[,TCR.cols[2]]))) |
      (is.na(TCR_data[,TCR.cols[3]]) & (!is.na(TCR_data[,TCR.cols[4]])))
  }

  return(TCR_data[,clone_definition_df$ColName, drop=F])
}

#' Build similarity matrices across all stringency levels
#' @param TCR_data Data frame with clone version columns.
#' @param ref_table Optional reference table to compare against.
#' @param clone_definition_df Clone definition data frame.
#' @return Named list of logical similarity matrices.
#' @keywords internal
#' @noRd
.get_similarity_matrices <- function(TCR_data, ref_table = NULL, clone_definition_df = .clone.definition.df){

  # Sanity check:
  .check_clone_definition_df(clone_definition_df)

  if (is.null(ref_table)){
    ref_table <- TCR_data
  } else {
    necessary.cols.ref <- c(clone_definition_df$ColName)
    if (!all(necessary.cols.ref %in% colnames(ref_table))){
      missing.cols <- necessary.cols.ref[!necessary.cols.ref %in% colnames(ref_table)]
      stop("Missing required columns in ref_table: ", paste(missing.cols, collapse = ", "),
           ". Expected: ", paste(necessary.cols.ref, collapse = ", "))
    }
  }

  # Get the empty matrices and fill them up!
  n_rows <- nrow(TCR_data)
  n_cols <- nrow(ref_table)
  mat.list <- list()
  for (clone_def in  unique(clone_definition_df$CloneDef)){
    mat.list[[clone_def]] <- matrix(data = 0, nrow = n_rows, ncol = n_cols)
  }

  column.select.list <- split(clone_definition_df, clone_definition_df$CloneDef)
  for (row_i in 1:n_rows){
    # Get ref clone:
    ref.clone <- TCR_data[row_i, ]$Sample_A1_A2_B1_B2
    # Get the matched dataframe
    match.df <- ref.clone == ref_table[, clone_definition_df$ColName]
    # Get the matched indexes
    for (clone_def in  unique(clone_definition_df$CloneDef)){
      mat.list[[clone_def]][row_i,] <- rowSums(match.df[,column.select.list[[clone_def]]$ColName], na.rm = TRUE) > 0
    }
  }

  # repair some mapping issues:
  for (mat.list.name in names(mat.list)){
    mat.list.name.ref <- unique(clone_definition_df[as.character(clone_definition_df$CloneDef) == mat.list.name,]$CloneDefRef)
    if (mat.list.name.ref %in% names(mat.list)){
      # Correct:
      mat.list[[mat.list.name]] <- mat.list[[mat.list.name]] + mat.list[[mat.list.name.ref]]
      mat.list[[mat.list.name]] <- mat.list[[mat.list.name]] > 0
      # Add:
      mat.list[[mat.list.name]] <- mat.list[[mat.list.name]] + mat.list[[mat.list.name.ref]]
    }
  }

  return(mat.list)
}

#' Derive clone IDs from adjacency matrix using graph components
#' @param TCR_data Data frame with TCR chain and Sample columns.
#' @param matrix_adj Adjacency matrix.
#' @param reference_cluster Optional reference cluster for conflict resolution.
#' @param group_by_sample_id Logical. If TRUE, clone IDs include sample prefix.
#' @return List with Clone.df (cluster and CloneID) and ref.table.
#' @keywords internal
#' @noRd
.get_cloneID_from_AdjMat <- function(TCR_data, matrix_adj, reference_cluster = NULL, group_by_sample_id = TRUE){

  # Get disconnected components from graph:
  g <- graph_from_adjacency_matrix(adjmatrix = matrix_adj, mode = "undirected", diag = FALSE)
  TCR_data$cluster <- components(g)$membership
  TCR_data$cluster[rowSums(matrix_adj) == 0] <- as.numeric(NA)

  # Deal with possible clone match conflicts 
  # Clone match conflicts are when a clone is defined by more than 2 TRA or TRB alleles
  check.issues <- TCR_data %>%
    filter(!is.na(.data$cluster)) %>%
    group_by(.data$cluster) %>%
    reframe(
      Sample = unique(.data$Sample),
      .unique_TCR_element(.data$TRA_1, .data$TRA_2, "TRA"),
      .unique_TCR_element(.data$TRB_1, .data$TRB_2, "TRB")
    ) %>%
    as.data.frame()
  check.issues$Error <- (check.issues$TRA_n_chain > 2) | (check.issues$TRB_n_chain > 2)

  if (any(check.issues$Error)){
    if (is.null(reference_cluster)){
      warning("[VDJigsaw] Conflicts detected in cluster assignment (>2 chains per type) but no reference cluster provided for resolution")
    } else{
      # Get the problematic clusters:
      mask.error <- TCR_data$cluster %in% subset(check.issues, check.issues$Error == TRUE)$cluster
      # Find the ones to replace:
      cluster.2.replace <- reference_cluster[mask.error]
      # Sanity check: Check if all clusters in the reference are defined in the mask
      if (length(cluster.2.replace) != sum(reference_cluster %in% unique(cluster.2.replace))){
        warning("[VDJigsaw] Conflicts detected and replacement clusters are not valid")
      }
      # Set the new clusters
      new.cluster.nb <- match(cluster.2.replace, unique(cluster.2.replace)) + max(check.issues$cluster, na.rm = TRUE)
      TCR_data[mask.error,]$cluster <- new.cluster.nb
    }
  }

  # Get clone ID
  if (group_by_sample_id){
    # To assign CloneID: Create a dense rank of clusters within each Sample
    orig.rownames <- rownames(TCR_data)
    TCR_data <- TCR_data %>%
      group_by(.data$Sample) %>%
      mutate(
        CloneID = if_else(
          condition = !is.na(.data$cluster),
          true = paste0("Clone_", .data$Sample, "_", dense_rank(.data$cluster)),
          false = as.character(NA)
          )
        ) %>%
      ungroup() %>%
      as.data.frame() %>%
      'rownames<-'(orig.rownames)
  } else {
    orig.rownames <- rownames(TCR_data)
    TCR_data <- TCR_data %>%
      mutate(
        CloneID = if_else(
          condition = !is.na(.data$cluster),
          true = paste0("Clone_", dense_rank(.data$cluster)),
          false = as.character(NA)
        )
      ) %>%
      as.data.frame() %>%
      'rownames<-'(orig.rownames)
  }

  # Get the reference table:
  ref_table <- TCR_data %>%
    filter(!is.na(.data$CloneID)) %>%
    group_by(.data$CloneID) %>%
    reframe(
      Sample = unique(.data$Sample),
      .unique_TCR_element(.data$TRA_1, .data$TRA_2, "TRA"),
      .unique_TCR_element(.data$TRB_1, .data$TRB_2, "TRB")
    ) %>%
    as.data.frame() %>%
    'rownames<-'(.$CloneID)

  # return results
  return.list <- list()
  return.list[["Clone.df"]] <- TCR_data[,c("cluster", "CloneID")]
  return.list[["ref.table"]] <- ref_table

  return(return.list)
}

#' Map cells to reference clonotypes using similarity matrix
#' @param similarity_mat Similarity matrix (cells x reference clonotypes).
#' @param ref_table Reference table with CloneID column.
#' @return Character vector of mapped clone IDs.
#' @keywords internal
#' @noRd
.get_cloneID_from_reference_and_AdjMat <- function(similarity_mat, ref_table){

  # Get the clone name from the matching elements
  map_res <- as.character(apply(X = similarity_mat, MARGIN = 1, FUN = function(x) ifelse(any(x > 0), paste0(ref_table$CloneID[x > 0], collapse = "|"), "")))
  map_res[map_res == ""] <- as.character(NA)

  # Identify conflicts:
  conflicting_maps <- rowSums(similarity_mat > 0) > 1
  if (any(conflicting_maps)){
    warning("[VDJigsaw] Found ", sum(conflicting_maps), " cells mapping to multiple clonotypes. ",
            "IDs are collapsed and separated by '|'")
  }

  return(map_res)
}

# ------------------------------------------------------------------------------
# identify_clonotypes
# ------------------------------------------------------------------------------

#' Identify clonotypes at multiple stringency levels
#'
#' Groups cells into clonotypes using a hierarchy of stringency levels, from
#' strict (all four chains must match) to permissive (a single allele suffices).
#' Uses graph-based clustering via igraph to find connected components of
#' matching cells. Supports parallel processing across samples.
#'
#' @param TCR_data A data frame with columns TRA_1, TRA_2, TRB_1, TRB_2
#'   containing validated TCR chain identifiers.
#' @param clone_definition_df Clone definition data frame specifying the
#'   stringency hierarchy. Defaults to the built-in definition.
#' @param num_cores Number of cores for parallel processing across samples.
#'   Set to \code{-1} to use all available cores. The actual number of cores
#'   used is capped at the number of available cores and the number of unique
#'   samples.
#' @param sample_col Name of the column to use as sample identifier. If
#'   \code{NULL}, all cells are treated as one sample.
#' @param clone_loose Which stringency level to use as the "loose" default.
#'   One of: "dual_chain_dual_allele", "dual_chain_one_partial",
#'   "dual_chain_both_partial", "single_chain_dual_allele",
#'   "single_chain_single_allele".
#' @param verbose Logical. If \code{TRUE}, print progress messages.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{TCR_data}{Data frame with clone ID columns for each stringency
#'       level (CloneID.dual_chain_dual_allele, etc.) plus CloneID.loose.}
#'     \item{ref_tables}{Named list of reference tables for each stringency
#'       level, each containing CloneID, Sample, and unique TCR chain info.}
#'   }
#'
#' @export
identify_clonotypes <- function(TCR_data,
                                clone_definition_df = .clone.definition.df,
                                num_cores = 1,
                                sample_col = NULL,
                                clone_loose = "single_chain_single_allele",
                                verbose = TRUE){

  # Part 1: Check validity of input data
  if (verbose) message("[VDJigsaw] Identifying clonotypes and building reference tables...")
  if (verbose) message("[VDJigsaw] Checking inputs...")
  .check_clone_definition_df(clone_definition_df)

  necessary.cols <- c("TRA_1", "TRA_2", "TRB_1", "TRB_2")
  if (!all(necessary.cols %in% colnames(TCR_data))){
    missing.cols <- necessary.cols[!necessary.cols %in% colnames(TCR_data)]
    stop("Missing required TCR chain columns in TCR_data: ", paste(missing.cols, collapse = ", "),
         ". Available columns: ", paste(colnames(TCR_data), collapse = ", "))
  }

  # Set up the sample id
  if (is.null(sample_col)){
    TCR_data$Sample <- "Sample"
  } else if (sample_col %in% colnames(TCR_data)){
    TCR_data$Sample <- TCR_data[[sample_col]]
  } else {
    stop("Sample column '", sample_col, "' not found in TCR_data. ",
         "Available columns: ", paste(colnames(TCR_data), collapse = ", "))
  }

  if (verbose) message("[VDJigsaw] Processing ", nrow(TCR_data), " cells across ",
                      length(unique(TCR_data$Sample)), " sample(s)")

  # Set up rownames to ensure correct mapping
  TCR_data$rownames.orig <- rownames(TCR_data)
  TCR_data$rownames.new <- sprintf(paste0("%0", (ceiling(log10(nrow(TCR_data)))+1), "d"), seq(1, nrow(TCR_data)))
  rownames(TCR_data) <- TCR_data$rownames.new

  # Set up parallel backend
  available_cores <- detectCores()
  n_samples <- length(unique(TCR_data$Sample))

  if (num_cores == -1) {
    num_cores <- available_cores
  }

  if (num_cores > available_cores) {
    warning("[VDJigsaw] Requested ", num_cores, " cores but only ", available_cores, " available. Using ", available_cores, " cores instead.")
    num_cores <- available_cores
  }

  num_cores <- min(num_cores, n_samples)

  if (verbose) message("[VDJigsaw] Using ", num_cores, " core(s) for parallel processing")

  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  # Export necessary functions and data to parallel workers
  # (This ensures the internal functions are available in each worker)
  clusterExport(cl,
                varlist = c(".get_clone_versions", ".get_similarity_matrices",
                           ".get_cloneID_from_AdjMat", ".unique_TCR_element",
                           ".check_clone_definition_df", ".clone.definition.df"),
                envir = environment())

  # Process each group in parallel
  grp <- NULL # to avoid R CMD check NOTE about undefined global variable
  results <- foreach(grp = unique(TCR_data$Sample), .packages = c("dplyr", "tidyr", "igraph")) %dopar% {
    # Process single group
    TCR_data_grp <- TCR_data[TCR_data$Sample == grp,]

    # Get clone versions
    TCR_data_grp <- cbind(TCR_data_grp, .get_clone_versions(TCR_data = TCR_data_grp))

    # Get similarity matrices
    similarity_matrices <- .get_similarity_matrices(TCR_data = TCR_data_grp)

    # Get CloneIDs
    ref_table_list <- list()
    for (mat.name in names(similarity_matrices)){
      matrix_adj <- similarity_matrices[[mat.name]]
      reference <- unique(clone_definition_df[as.character(clone_definition_df$CloneDef) == mat.name,]$CloneDefRef)
      if (reference %in% names(similarity_matrices)){
        reference_cluster <- TCR_data_grp[[paste0("cluster.", reference)]]
      } else {
        reference_cluster <- NULL
      }

      cloneID_res_list <- .get_cloneID_from_AdjMat(
        TCR_data = TCR_data_grp,
        matrix_adj = matrix_adj,
        reference_cluster = reference_cluster)
      ref_table_list[[mat.name]] <- cloneID_res_list$ref.table
      colnames(cloneID_res_list$Clone.df) <- paste0(colnames(cloneID_res_list$Clone.df), ".", mat.name)
      TCR_data_grp <- cbind(TCR_data_grp, cloneID_res_list$Clone.df[rownames(TCR_data_grp),])
    }

    # Return results for this group
    list(
      TCR_data = TCR_data_grp[, c(necessary.cols, paste0("CloneID.", names(similarity_matrices)))],
      ref_tables = ref_table_list
    )
  }

  stopCluster(cl)

  # Combine results from parallel processing
  if (verbose) message("[VDJigsaw] Combining results from parallel processing...")
  TCR_data_stacked <- data.frame()
  ref_tables <- list()
  for (i in seq_along(results)) {
    res_list_grp <- results[[i]]
    TCR_data_stacked <- rbind(TCR_data_stacked, res_list_grp$TCR_data)
    for (ref_table_name in names(res_list_grp$ref_tables)){
      if (ref_table_name %in% names(ref_tables)){
        ref_tables[[ref_table_name]] <- rbind(ref_tables[[ref_table_name]], res_list_grp$ref_tables[[ref_table_name]])
      } else {
        ref_tables[[ref_table_name]] <- res_list_grp$ref_tables[[ref_table_name]]
      }
    }
  }

  # Ensure rows are ordered as original TCR_data
  TCR_data_stacked <- TCR_data_stacked[rownames(TCR_data),]
  
  # Set back to old row.names
  rownames(TCR_data_stacked) <- TCR_data$rownames.orig

  # Add clone loose
  ref_tables[["loose"]] <- ref_tables[[clone_loose]]
  TCR_data_stacked$CloneID.loose <- TCR_data_stacked[[paste0("CloneID.", clone_loose)]]

  # Report summary
  if (verbose) {
    for (stringency in names(ref_tables)) {
      if (stringency != "loose") {
        n_clones <- nrow(ref_tables[[stringency]])
        message("[VDJigsaw]   ", stringency, ": ", n_clones, " clonotypes identified")
      }
    }
    message("[VDJigsaw] Clonotype identification complete!")
  }

  # Return results
  res_list <- list()
  res_list[["TCR_data"]] <- TCR_data_stacked
  res_list[["ref_tables"]] <- ref_tables

  return(res_list)
}

# ------------------------------------------------------------------------------
# map_clone_id
# ------------------------------------------------------------------------------

#' Map clone IDs from VDJigsaw results to wide-format data
#'
#' Extracts the clone ID columns from VDJigsaw results and aligns them to
#' the rows of the wide-format VDJ data frame.
#'
#' @param VDJigsaw_res Output from \code{\link{identify_clonotypes}}.
#' @param VDJ_contigs_wide Wide-format data frame (output of
#'   \code{\link{pivot_VDJ}}).
#' @param cols_to_map Character vector of clone ID column names to extract.
#'
#' @return A data frame with the requested clone ID columns, aligned to the
#'   rows of \code{VDJ_contigs_wide}.
#'
#' @export
map_clone_id <- function(VDJigsaw_res, VDJ_contigs_wide, cols_to_map = c(
  'CloneID.dual_chain_dual_allele', 'CloneID.loose', 'CloneID.dual_chain_one_partial', 'CloneID.dual_chain_both_partial',
  'CloneID.single_chain_dual_allele', 'CloneID.single_chain_single_allele')){
  VDJigsaw_res$TCR_data[rownames(VDJ_contigs_wide), cols_to_map]
}

# ------------------------------------------------------------------------------
# assign_clonotype (main pipeline)
# ------------------------------------------------------------------------------

#' Assign clonotypes from raw VDJ contig data
#'
#' Complete VDJigsaw pipeline that takes raw VDJ contig annotations, pivots
#' them to wide format, validates TCR chains, identifies clonotypes at multiple
#' stringency levels, and returns the annotated data with clone assignments.
#'
#' @param VDJ_data A data frame of VDJ contig annotations in 10x Genomics
#'   format.
#' @param sample_col Name of the column in \code{VDJ_data} to use as sample
#'   identifier. If \code{NULL}, all cells are assigned to "Sample".
#' @param clone_definition_df Clone definition data frame specifying the
#'   stringency hierarchy. Defaults to the built-in definition.
#' @param is_cell Logical. If \code{TRUE}, filter to rows where
#'   \code{is_cell == "true"}.
#' @param high_confidence Logical. If \code{TRUE}, filter to rows where
#'   \code{high_confidence == "true"}.
#' @param productive Logical. If \code{TRUE}, filter to rows where
#'   \code{productive == "true"}.
#' @param full_length Logical. If \code{TRUE}, filter to rows where
#'   \code{full_length == "true"}.
#' @param remove_invalid_VDJ Logical. If \code{TRUE}, set invalid VDJ
#'   sequences to \code{NA}.
#' @param cdr3_col Which CDR3 column to use for building the composite TCR
#'   chain identifier. Either \code{"cdr3_nt"} (nucleotide, default) or
#'   \code{"cdr3"} (amino acid).
#' @param num_cores Number of cores for parallel processing. Set to \code{-1}
#'   to use all available cores. Capped at the number of available cores and
#'   the number of unique samples.
#' @param clone_loose Which stringency level to use as the default "loose"
#'   clone definition.
#' @param verbose Logical. If \code{TRUE}, print progress messages.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{TCR_data}{Wide-format data frame with TCR chains, clone IDs at
#'       each stringency level, and a default CloneID column.}
#'     \item{ref_tables}{Named list of reference tables for each stringency
#'       level.}
#'   }
#'
#' @examples
#' \dontrun{
#' VDJ_contigs <- read.csv("filtered_contig_annotations.csv")
#' result <- assign_clonotype(VDJ_contigs, sample_col = "origin")
#' TCR_data <- result$TCR_data
#' ref_tables <- result$ref_tables
#' }
#'
#' @export
assign_clonotype <- function(
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
    verbose = TRUE){

  if (verbose) message("[VDJigsaw] ========================================")
  if (verbose) message("[VDJigsaw] Starting VDJigsaw clonotype assignment")
  if (verbose) message("[VDJigsaw] ========================================")

  # Pivot into wide format (rows = barcode)
  VDJ_contigs_wide <- pivot_VDJ(
    VDJ_data = VDJ_data,
    sample_col = sample_col,
    is_cell = is_cell,
    high_confidence = high_confidence,
    productive = productive,
    full_length = full_length,
    cdr3_col = cdr3_col,
    verbose = verbose)

  # Validate and correct problematic barcodes
  VDJ_contigs_wide <- validate_TCR(
    TCR_data = VDJ_contigs_wide,
    remove_invalid_VDJ = remove_invalid_VDJ,
    verbose = verbose)

  # Identify clonotypes
  VDJigsaw_res <- identify_clonotypes(
    TCR_data = VDJ_contigs_wide,
    clone_definition_df = clone_definition_df,
    num_cores = num_cores,
    sample_col = "Sample",
    clone_loose = clone_loose,
    verbose = verbose
  )

  # Align the clontype definition to the the wide table
  VDJigsaw_res_alligned <- map_clone_id(
    VDJigsaw_res = VDJigsaw_res,
    VDJ_contigs_wide = VDJ_contigs_wide)

  # Get the loose definition
  VDJigsaw_res_alligned$CloneID <- VDJigsaw_res_alligned[[paste0("CloneID.", "loose")]]

  # Map to the the wide table for complete results
  VDJ_contigs_wide <- cbind(VDJ_contigs_wide, VDJigsaw_res_alligned)

  # Re-order the columns
  column.order <- unique(c("Sample", "barcode", colnames(VDJigsaw_res_alligned), colnames(VDJ_contigs_wide)))
  VDJ_contigs_wide <- VDJ_contigs_wide[, column.order]

  # Prepare output (for simplicity I just modify `TCR_data` output of `VDJigsaw_res`)
  VDJigsaw_res$TCR_data <- VDJ_contigs_wide
  
  # Final messages
  if (verbose) {
    n_with_clones <- sum(!is.na(VDJ_contigs_wide$CloneID))
    n_total <- nrow(VDJ_contigs_wide)
    message("[VDJigsaw] ========================================")
    message("[VDJigsaw] Pipeline complete!")
    message("[VDJigsaw] Assigned ", n_with_clones, " of ", n_total, " cells (",
            round(100 * n_with_clones / n_total, 1), "%) to clonotypes")
    message("[VDJigsaw] ========================================")
  }

  return(VDJigsaw_res)
}
