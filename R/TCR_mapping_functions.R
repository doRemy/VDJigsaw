# TCR mapping functions
# @Author Rémy Pétremand

#' @importFrom dplyr %>% group_by ungroup mutate filter reframe summarise
#'   if_else dense_rank n select all_of across where group_by_at summarise
#'   .data mutate_at
#' @importFrom tidyr pivot_wider separate_wider_delim
#' @importFrom igraph graph_from_adjacency_matrix components
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores clusterExport
#' @importFrom stats setNames runif
#' @importFrom tibble tibble
#' @importFrom stringr regex
NULL

# Global variable declarations to avoid R CMD check NOTEs
utils::globalVariables(".")

# ------------------------------------------------------------------------------
# Default clone.definition.df and mapping levels
# ------------------------------------------------------------------------------

.clone.definition.df <- rbind(
  c("ColName",              "CloneDef",                  "CloneDefRef",               "TRA_1", "TRA_2", "TRB_1", "TRB_2"),
  # dual_chain_dual_allele combinations
  c("Sample_A1_A2_B1_B2", "dual_chain_dual_allele",    "None",                      "TRA_1", "TRA_2", "TRB_1", "TRB_2"),
  c("Sample_A1_A2_B2_B1", "dual_chain_dual_allele",    "None",                      "TRA_1", "TRA_2", "TRB_2", "TRB_1"),
  c("Sample_A2_A1_B1_B2", "dual_chain_dual_allele",    "None",                      "TRA_2", "TRA_1", "TRB_1", "TRB_2"),
  c("Sample_A2_A1_B2_B1", "dual_chain_dual_allele",    "None",                      "TRA_2", "TRA_1", "TRB_2", "TRB_1"),
  # dual_chain_one_partial combinations built from dual_chain_dual_allele
  c("Sample_A1_A2_B1_NA", "dual_chain_one_partial",    "dual_chain_dual_allele",    "TRA_1", "TRA_2", "TRB_1", "NAcol"),
  c("Sample_A1_A2_B2_NA", "dual_chain_one_partial",    "dual_chain_dual_allele",    "TRA_1", "TRA_2", "TRB_2", "NAcol"),
  c("Sample_A2_A1_B1_NA", "dual_chain_one_partial",    "dual_chain_dual_allele",    "TRA_2", "TRA_1", "TRB_1", "NAcol"),
  c("Sample_A2_A1_B2_NA", "dual_chain_one_partial",    "dual_chain_dual_allele",    "TRA_2", "TRA_1", "TRB_2", "NAcol"),
  c("Sample_A1_NA_B1_B2", "dual_chain_one_partial",    "dual_chain_dual_allele",    "TRA_1", "NAcol", "TRB_1", "TRB_2"),
  c("Sample_A1_NA_B2_B1", "dual_chain_one_partial",    "dual_chain_dual_allele",    "TRA_1", "NAcol", "TRB_2", "TRB_1"),
  c("Sample_A2_NA_B1_B2", "dual_chain_one_partial",    "dual_chain_dual_allele",    "TRA_2", "NAcol", "TRB_1", "TRB_2"),
  c("Sample_A2_NA_B2_B1", "dual_chain_one_partial",    "dual_chain_dual_allele",    "TRA_2", "NAcol", "TRB_2", "TRB_1"),
  # dual_chain_both_partial combinations built from dual_chain_one_partial
  c("Sample_A1_NA_B1_NA", "dual_chain_both_partial",    "dual_chain_one_partial",   "TRA_1", "NAcol", "TRB_1", "NAcol"),
  c("Sample_A1_NA_B2_NA", "dual_chain_both_partial",    "dual_chain_one_partial",   "TRA_1", "NAcol", "TRB_2", "NAcol"),
  c("Sample_A2_NA_B1_NA", "dual_chain_both_partial",    "dual_chain_one_partial",   "TRA_2", "NAcol", "TRB_1", "NAcol"),
  c("Sample_A2_NA_B2_NA", "dual_chain_both_partial",    "dual_chain_one_partial",   "TRA_2", "NAcol", "TRB_2", "NAcol"),
  # single_chain_dual_allele combinations built from dual_chain_both_partial
  c("Sample_A1_A2_NA_NA", "single_chain_dual_allele",   "dual_chain_both_partial",  "TRA_1", "TRA_2", "NAcol", "NAcol"),
  c("Sample_A2_A1_NA_NA", "single_chain_dual_allele",   "dual_chain_both_partial",  "TRA_2", "TRA_1", "NAcol", "NAcol"),
  c("Sample_NA_NA_B1_B2", "single_chain_dual_allele",   "dual_chain_both_partial",  "NAcol", "NAcol", "TRB_1", "TRB_2"),
  c("Sample_NA_NA_B2_B1", "single_chain_dual_allele",   "dual_chain_both_partial",  "NAcol", "NAcol", "TRB_2", "TRB_1"),
  # single_chain_single_allele combinations built from single_chain_dual_allele
  c("Sample_A1_NA_NA_NA", "single_chain_single_allele", "single_chain_dual_allele", "TRA_1", "NAcol", "NAcol", "NAcol"),
  c("Sample_A2_NA_NA_NA", "single_chain_single_allele", "single_chain_dual_allele", "TRA_2", "NAcol", "NAcol", "NAcol"),
  c("Sample_NA_NA_B1_NA", "single_chain_single_allele", "single_chain_dual_allele", "NAcol", "NAcol", "TRB_1", "NAcol"),
  c("Sample_NA_NA_B2_NA", "single_chain_single_allele", "single_chain_dual_allele", "NAcol", "NAcol", "TRB_2", "NAcol")
)
# Convert to a data frame for proper visualization
.clone.definition.df <- as.data.frame(.clone.definition.df, stringsAsFactors = FALSE)
colnames(.clone.definition.df) <- .clone.definition.df[1, ]
.clone.definition.df <- .clone.definition.df[-1, ]

# Set clone def as factor:
.mapping.levels <- c("dual_chain_dual_allele", "dual_chain_one_partial", "dual_chain_both_partial", "single_chain_dual_allele", "single_chain_single_allele")
.clone.definition.df$CloneDef <- factor(.clone.definition.df$CloneDef, levels = .mapping.levels)

# ------------------------------------------------------------------------------
# internal helpers
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

# ------------------------------------------------------------------------------
# Get data
# ------------------------------------------------------------------------------

#' Pivot VDJ contig data to wide format
#'
#' Converts VDJ contig annotation data from long format (one row per contig)
#' to wide format (one row per cell barcode), creating composite TCR chain
#' identifiers (TRA_1, TRA_2, TRB_1, TRB_2) from V gene, CDR3, and J gene.
#'
#' @param VDJ_data A data frame of VDJ contig annotations in 10x Genomics
#'   format. Must contain columns: barcode, chain, v_gene, cdr3, j_gene, umis,
#'   and other standard 10x VDJ columns.
#' @param sample_col Name of the column in \code{VDJ_data} to use as sample
#'   identifier. If \code{NULL}, all cells are assigned to a single sample
#'   called "Sample".
#' @param is_cell Logical. If \code{TRUE}, filter to rows where
#'   \code{is_cell == "true"}.
#' @param high_confidence Logical. If \code{TRUE}, filter to rows where
#'   \code{high_confidence == "true"}.
#' @param productive Logical. If \code{TRUE}, filter to rows where
#'   \code{productive == "true"}.
#' @param full_length Logical. If \code{TRUE}, filter to rows where
#'   \code{full_length == "true"}.
#' @param verbose Logical. If \code{TRUE}, print progress messages.
#'
#' @return A data frame in wide format with one row per cell barcode. Contains
#'   columns for sample, barcode, TCR chain identifiers (TRA_1, TRA_2, TRB_1,
#'   TRB_2), contig metadata, and individual chain component columns.
#'
#' @export
pivot_VDJ <- function(VDJ_data, sample_col = NULL, is_cell = FALSE, high_confidence = FALSE, productive = FALSE, full_length = FALSE, verbose = TRUE){

  if (verbose) message("[VDJigsaw] Converting VDJ data to wide format...")
  n_original <- nrow(VDJ_data)

  # Filtering the barcode
  if (is_cell){
    VDJ_data <- VDJ_data[VDJ_data$is_cell == "true",]
    if (verbose) message("[VDJigsaw]   Filtered to cells marked as is_cell: ", nrow(VDJ_data), " rows")
  }
  if (high_confidence){
    VDJ_data <- VDJ_data[VDJ_data$high_confidence == "true",]
    if (verbose) message("[VDJigsaw]   Filtered to high confidence: ", nrow(VDJ_data), " rows")
  }
  if (productive){
    VDJ_data <- VDJ_data[VDJ_data$productive == "true",]
    if (verbose) message("[VDJigsaw]   Filtered to productive: ", nrow(VDJ_data), " rows")
  }
  if (full_length){
    VDJ_data <- VDJ_data[VDJ_data$full_length == "true",]
    if (verbose) message("[VDJigsaw]   Filtered to full length: ", nrow(VDJ_data), " rows")
  }

  # get the Sample column
  if (is.null(sample_col)){
    VDJ_data$Sample <- "Sample"
  } else if (sample_col %in% colnames(VDJ_data)){
    VDJ_data$Sample <- VDJ_data[[sample_col]]
  } else {
    stop("Sample column '", sample_col, "' not found in VDJ_data. Available columns: ", paste(colnames(VDJ_data), collapse = ", "))
  }


  # Keep the columns of interest
  cell.info.col <- c('Sample', 'barcode')
  chain.info.col <- c(
    'contig_id', 'length',
    'raw_clonotype_id', 'raw_consensus_id', 'exact_subclonotype_id',
    'chain', 'v_gene', 'd_gene', 'j_gene', 'c_gene',
    'reads', 'umis',
    'fwr1', 'fwr1_nt', 'cdr1', 'cdr1_nt',
    'fwr2', 'fwr2_nt', 'cdr2', 'cdr2_nt',
    'fwr3', 'fwr3_nt', 'cdr3', 'cdr3_nt',
    'fwr4', 'fwr4_nt'
  )
  VDJ_data <- VDJ_data[, c(cell.info.col, chain.info.col)]

  # Get the VDJ data wide format (row = barcode)
  if (verbose) message("[VDJigsaw]   Pivoting to wide format...")
  VDJ_data_wide <- VDJ_data %>%
    group_by(.data$barcode, .data$chain) %>%
    mutate(rank_alleles = rank(-.data$umis, ties.method = "first")) %>%
    ungroup() %>%
    pivot_wider(
      id_cols = all_of(cell.info.col),
      names_from = all_of(c("chain", "rank_alleles")),
      values_from = all_of(chain.info.col)) %>%
    as.data.frame()

  # sanity check
  if (any(duplicated(VDJ_data_wide$barcode))){
    warning("[VDJigsaw] Warning: Duplicated barcodes found after pivoting")
  }

  if (verbose) message("[VDJigsaw]   Wide format: ", nrow(VDJ_data_wide), " unique barcodes")

  # Get "TRA_1", "TRA_2", "TRB_1", "TRB_2" columns
  # Get the definition of "TRA_1", "TRA_2", "TRB_1", "TRB_2" columns
  if (verbose) message("[VDJigsaw]   Creating TRA/TRB chain columns...")
  process_col_list <- c("TRA_1", "TRA_2", "TRB_1", "TRB_2")
  names(process_col_list) <- process_col_list
  process_col_list <- lapply(
    X = process_col_list,
    FUN = function(x) paste(c("v_gene", "cdr3", "j_gene"), x, sep = "_"))
  # Processing
  for (col_id in names(process_col_list)){
    VDJ_data_wide[[col_id]] <- apply(
      X = VDJ_data_wide[, process_col_list[[col_id]]],
      MARGIN = 1,
      FUN = function(x) paste(x, sep = "_", collapse = "_"))
    mask.invalid <- apply(
      X = VDJ_data_wide[, process_col_list[[col_id]]],
      MARGIN = 1,
      FUN = function(x) all(is.na(x)))
    VDJ_data_wide[[col_id]][mask.invalid] <- as.character(NA)
  }

  # Get "contig_id" and clonotype columns
  # Get the definition of 'contig_id', 'raw_clonotype_id', 'raw_consensus_id', 'exact_subclonotype_id'
  process_col_list <- list('contig_id', 'raw_clonotype_id', 'raw_consensus_id', 'exact_subclonotype_id')
  names(process_col_list) <- process_col_list
  process_col_list <- lapply(
    X = process_col_list,
    FUN = function(x) paste(x, c("TRA_1", "TRA_2", "TRB_1", "TRB_2"), sep = "_"))
  # Processing
  for (col_id in names(process_col_list)){
    VDJ_data_wide[[paste(col_id, "combined", sep = "_")]] <- apply(
      X = VDJ_data_wide[, process_col_list[[col_id]]],
      MARGIN = 1,
      FUN = function(x) paste(unique(x[!is.na(x)]), sep = " ", collapse = " "))
    VDJ_data_wide[[paste(col_id, "n", sep = "_")]] <- apply(
      X = VDJ_data_wide[, process_col_list[[col_id]]],
      MARGIN = 1,
      FUN = function(x) length(unique(x[!is.na(x)])))
  }

  # Re.order data and keep columns of interst
  column.order <- c(
    'Sample', 'barcode',
    'TRA_1', 'TRA_2', 'TRB_1', 'TRB_2',
    'contig_id_combined', 'contig_id_n',
    'raw_clonotype_id_combined', 'raw_clonotype_id_n',
    'raw_consensus_id_combined', 'raw_consensus_id_n',
    'exact_subclonotype_id_combined', 'exact_subclonotype_id_n',
    unlist(lapply(
      X = c(
        'length',
        'v_gene', 'd_gene', 'j_gene', 'c_gene',
        'reads', 'umis',
        'fwr1', 'fwr1_nt', 'cdr1', 'cdr1_nt',
        'fwr2', 'fwr2_nt', 'cdr2', 'cdr2_nt',
        'fwr3', 'fwr3_nt', 'cdr3', 'cdr3_nt',
        'fwr4', 'fwr4_nt'),
      FUN = function(x) paste(x, c("TRA_1", "TRA_2", "TRB_1", "TRB_2"), sep = "_")))
  )
  VDJ_data_wide <- VDJ_data_wide[, column.order]

  if (verbose) message("[VDJigsaw] Wide format conversion complete: ", ncol(VDJ_data_wide), " columns")

  return(VDJ_data_wide)
}

#' Validate and correct TCR chain data
#'
#' Validates TCR chain columns (TRA_1, TRA_2, TRB_1, TRB_2) by checking
#' V-CDR3-J format, reordering chains into the canonical V_CDR3_J format,
#' removing duplicate alleles (where allele 1 equals allele 2), and swapping
#' alleles when allele 2 is defined but allele 1 is missing.
#'
#' @param TCR_data A data frame with columns TRA_1, TRA_2, TRB_1, TRB_2
#'   containing TCR chain identifiers in V_CDR3_J format.
#' @param remove_invalid_VDJ Logical. If \code{TRUE}, set invalid VDJ sequences
#'   (those with missing V, CDR3, or J components) to \code{NA}. If
#'   \code{FALSE}, keep them as-is and flag them.
#' @param verbose Logical. If \code{TRUE}, print progress messages.
#'
#' @return The input data frame with corrected TCR chain columns and additional
#'   logical columns (TRA_1.invalid, TRA_2.invalid, TRB_1.invalid,
#'   TRB_2.invalid) indicating which chains had invalid VDJ format.
#'
#' @export
validate_TCR <- function(TCR_data, remove_invalid_VDJ = FALSE, verbose = TRUE){

  if (verbose) message("[VDJigsaw] Validating and correcting TCR data...")

  necessary.cols <- c("TRA_1", "TRA_2", "TRB_1", "TRB_2")
  if (!all(necessary.cols %in% colnames(TCR_data))){
    missing.cols <- necessary.cols[!necessary.cols %in% colnames(TCR_data)]
    stop("Missing required columns in TCR_data: ", paste(missing.cols, collapse = ", "),
         ". Available columns: ", paste(colnames(TCR_data), collapse = ", "))
  }

  # Check that every chain is in the v-d-j format and that they do not have missing values
  n_invalid_total <- 0
  for (TCR_chain in c("TRA_1", "TRA_2", "TRB_1", "TRB_2")){
    # Part 1: Get chain parts
    TR_VDJ_chain <- TCR_data[, TCR_chain, drop = F] %>%
      separate_wider_delim(cols = all_of(TCR_chain), delim = "_",  names = c("A", "B", "C"), cols_remove = FALSE) %>%
      mutate(across(where(is.character), trimws)) %>%
      mutate(across(where(is.character), function(x){ifelse(x == "", as.character(NA), x)})) %>%
      as.data.frame()

    # Part 2: Identify which part is which
    col.TRV <- names(which.max(apply(X = TR_VDJ_chain[!is.na(TR_VDJ_chain[[TCR_chain]]), c("A", "B", "C")], MARGIN = 2, FUN = function(x){sum(grepl("TR[AB]V", x))}, simplify = F)))
    col.TRJ <- names(which.max(apply(X = TR_VDJ_chain[!is.na(TR_VDJ_chain[[TCR_chain]]), c("A", "B", "C")], MARGIN = 2, FUN = function(x){sum(grepl("TR[AB]J", x))}, simplify = F)))
    col.CDR <- c("A", "B", "C")[!(c("A", "B", "C") %in% c(col.TRV, col.TRJ))]
    TR_VDJ_chain$TRV <- TR_VDJ_chain[[col.TRV]]
    TR_VDJ_chain$TRV[TR_VDJ_chain$TRV == "NA"] <- NA
    TR_VDJ_chain$TRJ <- TR_VDJ_chain[[col.TRJ]]
    TR_VDJ_chain$TRJ[TR_VDJ_chain$TRJ == "NA"] <- NA
    TR_VDJ_chain$CDR <- TR_VDJ_chain[[col.CDR]]
    TR_VDJ_chain$CDR[TR_VDJ_chain$CDR == "NA"] <- NA

    # Part 3: Get VDJ
    TR_VDJ_chain$VDJ <- paste(TR_VDJ_chain$TRV, TR_VDJ_chain$CDR, TR_VDJ_chain$TRJ, sep = "_")
    TR_VDJ_chain$VDJ[is.na(TR_VDJ_chain[[TCR_chain]])] <- as.character(NA)

    # Part 4: Identify invalid VDJ sequences i.e. VDJ sequences with missing part


    TR_VDJ_chain$Invalid.VDJ <- !is.na(TR_VDJ_chain$VDJ) & ((rowSums(is.na(TR_VDJ_chain[, c("TRV", "CDR", "TRJ")])) %% 3) != 0)
    if (any(TR_VDJ_chain$Invalid.VDJ)){
      n_invalid <- sum(TR_VDJ_chain$Invalid.VDJ)
      n_invalid_total <- n_invalid_total + n_invalid
      warning.message <- paste0("[VDJigsaw] Found ", n_invalid, " invalid ", TCR_chain, " chains with missing V-D-J parts ",
                                "(see ", TCR_chain, ".invalid column)")
      if (remove_invalid_VDJ){
        if (verbose) message(warning.message, " - setting to NA")
        TR_VDJ_chain$VDJ[TR_VDJ_chain$Invalid.VDJ] <- as.character(NA)
      } else {
        if (verbose) message(warning.message, " - keeping as-is (set remove_invalid_VDJ=TRUE to remove)")
      }
    }

    # Final part: Keep important data
    TCR_data[[TCR_chain]] <- TR_VDJ_chain$VDJ
    TCR_data[[paste0(TCR_chain, ".invalid")]] <- TR_VDJ_chain$Invalid.VDJ
  }

  if (verbose && n_invalid_total > 0) {
    message("[VDJigsaw] Total invalid chains across all chain types: ", n_invalid_total)
  }

  # Check that the second allele is NA if both are defined
  n_corrections <- 0
  mask.equal <- TCR_data$TRA_1 == TCR_data$TRA_2; mask.equal[is.na(mask.equal)] <- FALSE
  if (any(mask.equal)){
    n_corrections <- n_corrections + sum(mask.equal)
    TCR_data[mask.equal,]$TRA_2 <- as.character(NA)
  }
  mask.equal <- TCR_data$TRB_1 == TCR_data$TRB_2; mask.equal[is.na(mask.equal)] <- FALSE
  if (any(mask.equal)){
    n_corrections <- n_corrections + sum(mask.equal)
    TCR_data[mask.equal,]$TRB_2 <- as.character(NA)
  }

  if (verbose && n_corrections > 0) {
    message("[VDJigsaw] Removed ", n_corrections, " duplicate alleles (allele 1 == allele 2)")
  }

  # If allele 2 is defined and not 1...
  n_swaps <- 0
  mask.1 <- is.na(TCR_data$TRA_1) & !is.na(TCR_data$TRA_2)
  if (any(mask.1)){
    n_swaps <- n_swaps + sum(mask.1)
    TCR_data$TRA_1[mask.1] <- TCR_data$TRA_2[mask.1]
    TCR_data$TRA_2[mask.1] <- as.character(NA)
  }

  mask.1 <- is.na(TCR_data$TRB_1) & !is.na(TCR_data$TRB_2)
  if (any(mask.1)){
    n_swaps <- n_swaps + sum(mask.1)
    TCR_data$TRB_1[mask.1] <- TCR_data$TRB_2[mask.1]
    TCR_data$TRB_2[mask.1] <- as.character(NA)
  }

  if (verbose && n_swaps > 0) {
    message("[VDJigsaw] Swapped ", n_swaps, " alleles where allele_2 was defined but allele_1 was not")
  }

  if (verbose) message("[VDJigsaw] TCR data correction complete")

  return(TCR_data)
}

# ------------------------------------------------------------------------------
# VDJigsaw internal functions and methods
# ------------------------------------------------------------------------------

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
    # Set the name of the new column
    new.col <- clone_definition_df[row_i, "ColName"]
    # Set the TCR columns
    TCR.cols <- as.character(clone_definition_df[row_i, c("TRA_1", "TRA_2", "TRB_1", "TRB_2")])
    # Set the clone column definition
    clone.cols <- c("Sample", TCR.cols)

    # Get the clone version
    TCR_data[[new.col]] <- apply(X = TCR_data[,clone.cols], MARGIN = 1, FUN = paste, sep = "_", collapse = "_")

    # Remove the aberrant ones:
    # <Sample>_NA_NA_NA_NA to NA
    TCR_data[[new.col]][rowSums(is.na(TCR_data[,c("TRA_1", "TRA_2", "TRB_1", "TRB_2")])) == 4] <- as.character(NA)

    # When the second allele is defined and the first one not --> Aberrant
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

  # Deal with possible conficts i.e. when one TCR match with antoher that it should not.
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
    # Create a dense rank of clusters within each Sample to get the Clone IDs
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
  if (num_cores >= detectCores()){
    num_cores <- detectCores() - 1
    warning("[VDJigsaw] Requested ", num_cores, " cores but only ", detectCores(),
            " available. Using ", num_cores, " cores instead")
  }

  if (verbose) message("[VDJigsaw] Using ", num_cores, " core(s) for parallel processing")

  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  # Export necessary functions and data to parallel workers
  # This ensures the internal functions are available in each worker
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

  # Return parallel results
  res_list <- list()
  res_list[["TCR_data"]] <- TCR_data_stacked
  res_list[["ref_tables"]] <- ref_tables

  return(res_list)
}

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
# VDJigsaw pipeline
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
#' @param num_cores Number of cores for parallel processing.
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

  # Map to the wide table
  VDJ_contigs_wide <- map_clone_id(
    VDJigsaw_res = VDJigsaw_res,
    VDJ_contigs_wide = VDJ_contigs_wide)

  # Get the loose definition
  VDJ_contigs_wide$CloneID <- VDJ_contigs_wide[[paste0("CloneID.", "loose")]]

  # Prepare output
  VDJigsaw_res$TCR_data <- VDJ_contigs_wide

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

# ------------------------------------------------------------------------------
# Visualization
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
  TCR_data <- TCR_data %>% group_by_at(clone.columns) %>% summarise(n = n(), .groups = "keep") %>% as.data.frame()
  TCR_data[,track.change.columns] <- matrix(FALSE, nrow = nrow(TCR_data), ncol = length(track.change.columns))
  TCR_data[,track.change.columns[1]] <- TRUE

  # Part 3: Track change
  for (i in 2:length(clone.columns)){
    count.df <- TCR_data %>%
      group_by(.data[[clone.columns[i]]]) %>%
      summarise(n_dup = length(unique(.data[[clone.columns[i-1]]])), .groups = "keep") %>%
      filter(.data$n_dup > 1) %>%
      as.data.frame()
    mask.dup <- TCR_data[[clone.columns[i]]] %in% count.df[[clone.columns[i]]]
    TCR_data[mask.dup, track.change.columns[i]] <- TRUE
  }

  # Part 4: Get the summary and return it
  summary.df <- TCR_data %>% group_by_at(track.change.columns) %>% summarise(n_clones = n(), .groups = "keep") %>% as.data.frame()

  return(summary.df)
}

# ------------------------------------------------------------------------------
# Mapping functions
# ------------------------------------------------------------------------------

#' Assign clonotypes by mapping to a reference table
#'
#' Maps TCR data to an existing reference clonotype table using similarity
#' matrices at a specified stringency level. Useful for mapping new data
#' against previously identified clonotypes.
#'
#' @param TCR_data A data frame with columns Sample, TRA_1, TRA_2, TRB_1,
#'   TRB_2.
#' @param ref_table A reference table with columns Sample, TRA_1, TRA_2,
#'   TRB_1, TRB_2, and CloneID (as output by \code{\link{identify_clonotypes}}).
#' @param similarity_matrix_type Which stringency level to use for matching.
#'   One of the values in \code{.mapping.levels}.
#' @param clone_definition_df Clone definition data frame.
#'
#' @return A character vector of mapped clone IDs, with \code{NA} for
#'   unmapped cells and pipe-separated IDs for cells matching multiple
#'   clonotypes.
#'
#' @export
assign_clonotype_from_reference <- function(TCR_data, ref_table, similarity_matrix_type = .mapping.levels, clone_definition_df = .clone.definition.df){

  # Part 1: Check validity of input data
  .check_clone_definition_df(clone_definition_df)

  similarity_matrix_type <- match.arg(similarity_matrix_type)

  necessary.cols <- c("Sample", "TRA_1", "TRA_2", "TRB_1", "TRB_2")
  if (!all(necessary.cols %in% colnames(TCR_data))){
    missing.cols <- necessary.cols[!necessary.cols %in% colnames(TCR_data)]
    stop("Missing required columns in TCR_data: ", paste(missing.cols, collapse = ", "),
         ". Available columns: ", paste(colnames(TCR_data), collapse = ", "))
  }
  if (!all(necessary.cols %in% colnames(ref_table))){
    missing.cols <- necessary.cols[!necessary.cols %in% colnames(ref_table)]
    stop("Missing required columns in ref_table: ", paste(missing.cols, collapse = ", "),
         ". Available columns: ", paste(colnames(ref_table), collapse = ", "))
  }

  # Part 2: Get the clone's possible versions
  TCR_data <- cbind(TCR_data, .get_clone_versions(TCR_data = TCR_data))
  ref_table <- cbind(ref_table, .get_clone_versions(TCR_data = ref_table))

  # Part 3: Get similarity matrices from the TCR data
  similarity_matrices <- .get_similarity_matrices(TCR_data = TCR_data, ref_table = ref_table)
  similarity_mat <- similarity_matrices[[similarity_matrix_type]]
  if (ncol(similarity_mat) != nrow(ref_table)){
    stop("Dimension mismatch in get_clone_id_from_reference(): ",
         "similarity matrix has ", ncol(similarity_mat), " columns but ref_table has ", nrow(ref_table), " rows")
  }

  # Part 4: get the mapping:
  TCR_data$CloneID.map <- .get_cloneID_from_reference_and_AdjMat(similarity_mat = similarity_mat, ref_table = ref_table)

  return(TCR_data$CloneID.map)
}


#' Annotate clone IDs with reference table metadata
#'
#' Given a vector of clone IDs (possibly pipe-separated for multi-mapped
#' cells), retrieves corresponding annotations from a reference table.
#' Handles numeric, logical, factor, and character columns appropriately.
#'
#' @param clone_ids Character vector of clone IDs. Pipe-separated IDs
#'   indicate cells mapping to multiple clonotypes.
#' @param ref_table A reference table with clone IDs as row names and
#'   annotation columns.
#'
#' @return A data frame with the same columns as \code{ref_table}, one row
#'   per input clone ID.
#'
#' @export
annotate_by_clone <- function(clone_ids, ref_table){

  # Part 1: Clean and Sanity checks
  TCR_data <- data.frame(CloneIDs = clone_ids)

  # Part 2: Separate the clone.col into separate clone definition
  max_nb_clones <- max(sapply(TCR_data$CloneIDs, function(x){length(strsplit(x = x, split = "[|]")[[1]])}))
  clone.sep.names <- paste0("C", seq(1, max_nb_clones))
  if (max_nb_clones == 1){
    TCR_data_sep <- TCR_data %>%
    mutate(C1 = .data$CloneIDs) %>%
    select(all_of(clone.sep.names)) %>%
    as.data.frame()
  } else {
    TCR_data_sep <- TCR_data %>%
      separate_wider_delim(
        cols = c("CloneIDs"),
        delim = regex("[|]"),
        names = clone.sep.names,
        too_few = "align_start",
        cols_remove = FALSE) %>%
      select(all_of(clone.sep.names)) %>%
      as.data.frame()
  }

  # Part 3: Mapping and obtaining one value per map:
  # Get the rows to remove later:
  mask.NA <- TCR_data_sep %>% mutate_at(clone.sep.names, function(x) x %in% rownames(ref_table))
  mask.NA <- !apply(mask.NA, 1, any)

  for (col_ in colnames(ref_table)){
    value.mapped <- TCR_data_sep %>% mutate_at(clone.sep.names, function(x){ref_table[x, col_]})

    if (is.numeric(ref_table[[col_]])){
      value.mapped <- apply(value.mapped, 1, function(x){ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE))})
    } else if (is.logical(ref_table[[col_]])){
      value.mapped <- apply(value.mapped, 1, any, na.rm = TRUE)
    } else if (is.factor(ref_table[[col_]])){
      values.ordered <- levels(ref_table[[col_]])
      value.mapped <- apply(value.mapped, 1, function(x){ifelse(any(values.ordered %in% x), values.ordered[values.ordered %in% x][1], NA)})
    } else {
      value.mapped <- apply(value.mapped, 1, function(x){unique(x)[1]})
    }
    value.mapped[mask.NA] <- NA

    TCR_data[[col_]] <- value.mapped
  }

  # Part 4: return the important columns:
  return(TCR_data[,colnames(ref_table)])
}

#' Map clonotypes from allele format to paired TCR format
#'
#' Maps TCR data from individual allele format (TRA_1/TRA_2/TRB_1/TRB_2) to
#' paired alpha/beta chain format (TRA/TRB). Useful when the reference table
#' uses only TRA and TRB (no second allele) and you want to map allele-level
#' data to this paired format.
#'
#' @param data_to_map A data frame with columns Sample, TRA_1, TRA_2, TRB_1,
#'   TRB_2.
#' @param ref_table A reference table with columns Sample, TRA, TRB, and the
#'   column specified by \code{col_to_map}.
#' @param col_to_map Name of the column in \code{ref_table} to map values
#'   from.
#' @param map_to_label Optional data frame for additional label mapping. Row
#'   names should correspond to mapped values. Factor and non-factor columns
#'   are handled differently.
#'
#' @return A data frame with mapping results including strict and loose
#'   matches, mapping type (strict/loose_A/loose_B/ND), and optionally
#'   label-mapped columns.
#'
#' @export
map_clonotypes_to_paired_TCR <- function(data_to_map, ref_table, col_to_map, map_to_label = NULL){

  # Define some functions:
  .func_map_val <- function(x){
    x = unique(x[x != "ND"])
    return(paste(x, collapse = "|"))
  }
  .func_map_val_to_label_factor <- function(x, map_to_label, col_label){
    vals.levels <- levels(map_to_label[[col_label]])
    vals.present <- c(vals.levels[vals.levels %in% map_to_label[x, col_label]])
    res <- paste(vals.present, collapse = "|")
    return(res)
  }
  .func_map_val_to_label <- function(x, map_to_label, col_label){
    vals.present <- unique(as.character(map_to_label[x, col_label]))
    vals.present[vals.present %in% c("NA", "ND")] <- as.character(NA)
    vals.present <- vals.present[!is.na(vals.present)]
    res <- paste(vals.present, collapse = "|")
    if (res == ""){res <- as.character(NA)}
    return(res)
  }

  # Check columns of data_to_map
  necessary.cols <- c("Sample", "TRA_1", "TRA_2", "TRB_1", "TRB_2")
  if (!all(necessary.cols %in% colnames(data_to_map))){
    missing.cols <- necessary.cols[!necessary.cols %in% colnames(data_to_map)]
    stop("Missing required columns in data_to_map: ", paste(missing.cols, collapse = ", "),
         ". Available columns: ", paste(colnames(data_to_map), collapse = ", "))
  }
  # Check columns of ref_table
  necessary.cols <- c("Sample", "TRA", "TRB", col_to_map)
  if (!all(necessary.cols %in% colnames(ref_table))){
    missing.cols <- necessary.cols[!necessary.cols %in% colnames(ref_table)]
    stop("Missing required columns in ref_table: ", paste(missing.cols, collapse = ", "),
         ". Available columns: ", paste(colnames(ref_table), collapse = ", "))
  }

  # Set the value_to_map
  ref_table$value_to_map <- ref_table[[col_to_map]]

  # Part 1: Pre-process
  data_to_map$Sample_TRA1_TRB1 <- paste(data_to_map$Sample, data_to_map$TRA_1, data_to_map$TRB_1, sep = "_")
  data_to_map$Sample_TRA1_TRB2 <- paste(data_to_map$Sample, data_to_map$TRA_1, data_to_map$TRB_2, sep = "_")
  data_to_map$Sample_TRA2_TRB1 <- paste(data_to_map$Sample, data_to_map$TRA_2, data_to_map$TRB_1, sep = "_")
  data_to_map$Sample_TRA2_TRB2 <- paste(data_to_map$Sample, data_to_map$TRA_2, data_to_map$TRB_2, sep = "_")
  data_to_map$Sample_TRA1 <-      paste(data_to_map$Sample, data_to_map$TRA_1, sep = "_")
  data_to_map$Sample_TRA2 <-      paste(data_to_map$Sample, data_to_map$TRA_2, sep = "_")
  data_to_map$Sample_TRB1 <-      paste(data_to_map$Sample, data_to_map$TRB_1, sep = "_")
  data_to_map$Sample_TRB2 <-      paste(data_to_map$Sample, data_to_map$TRB_2, sep = "_")

  ref_table$Sample_TRA_TRB <- paste(ref_table$Sample, ref_table$TRA, ref_table$TRB, sep = "_")
  ref_table$Sample_TRA <-     paste(ref_table$Sample, ref_table$TRA, sep = "_")
  ref_table$Sample_TRB <-     paste(ref_table$Sample, ref_table$TRB, sep = "_")

  # Part 2: Get mapping tables
  ref_TRA <- ref_table %>%
    filter(!is.na(.data$Sample_TRA)) %>%
    group_by(.data$Sample_TRA) %>%
    summarise(
      n_clones = n(),
      value_to_map = .func_map_val(.data$value_to_map)
    ) %>%
    as.data.frame() %>%
    'rownames<-'(.$Sample_TRA)

  ref_TRB <- ref_table %>%
    filter(!is.na(.data$Sample_TRB)) %>%
    group_by(.data$Sample_TRB) %>%
    summarise(
      n_clones = n(),
      value_to_map = .func_map_val(.data$value_to_map)
    ) %>%
    as.data.frame() %>%
    'rownames<-'(.$Sample_TRB)

  rownames(ref_table) <- ref_table$Sample_TRA_TRB

  # Part 3: Mapping
  data_map <- cbind(
    map_A1_B1 = ref_table[data_to_map$Sample_TRA1_TRB1, ]$value_to_map,
    map_A1_B2 = ref_table[data_to_map$Sample_TRA1_TRB2, ]$value_to_map,
    map_A2_B1 = ref_table[data_to_map$Sample_TRA2_TRB1, ]$value_to_map,
    map_A2_B2 = ref_table[data_to_map$Sample_TRA2_TRB2, ]$value_to_map,
    map_A1 = ref_TRA[data_to_map$Sample_TRA1,]$value_to_map,
    map_A2 = ref_TRA[data_to_map$Sample_TRA2,]$value_to_map,
    map_B1 = ref_TRB[data_to_map$Sample_TRB1,]$value_to_map,
    map_B2 = ref_TRB[data_to_map$Sample_TRB2,]$value_to_map,
    TRA_isna = is.na(data_to_map$TRA_1),
    TRB_isna = is.na(data_to_map$TRB_1)
  ) %>%
    as.data.frame()  %>%
    mutate(
      TRA_isna = as.logical(.data$TRA_isna),
      TRB_isna = as.logical(.data$TRB_isna)
    ) %>%
    'rownames<-'(rownames(data_to_map))
  data_map[is.na(data_map)] <- "ND"

  # Map the clone in a strict manner:
  map.AB.cols <- c("map_A1_B1", "map_A1_B2", "map_A2_B1", "map_A2_B2")
  data_map$map.strict <- apply(X = data_map[, map.AB.cols], MARGIN = 1, FUN = .func_map_val)

  # Map the clone in a loose manner:
  map.A.cols <- c("map_A1", "map_A2")
  data_map$map.A <- apply(X = data_map[, map.A.cols], MARGIN = 1, FUN = .func_map_val)

  map.B.cols <- c("map_B1", "map_B2")
  data_map$map.B <- apply(X = data_map[, map.B.cols], MARGIN = 1, FUN = .func_map_val)

  data_map$map.loose <- data_map$map.strict
  mask.map.strict.missing <- data_map$map.strict == ""
  mask.map.A.present <- data_map$map.A != ""
  mask.map.B.present <- data_map$map.B != ""

  mask.fill.mapA <- mask.map.strict.missing & mask.map.A.present & data_map$TRB_isna
  mask.fill.mapB <- mask.map.strict.missing & mask.map.B.present & data_map$TRA_isna

  # Note: TCR-beta chain has priority in loose matching
  data_map$map.loose[mask.fill.mapA] <- data_map$map.A[mask.fill.mapA]
  data_map$map.loose[mask.fill.mapB] <- data_map$map.B[mask.fill.mapB]


  # Part 4: Get some summary
  data_map$map.type <- "ND"
  data_map$map.type[!mask.map.strict.missing] <- "strict"
  data_map$map.type[mask.fill.mapA] <- "loose_A"
  data_map$map.type[mask.fill.mapB] <- "loose_B"

  # Part 5: Continue the mapping if the map_to_label is not NULL

  if (!is.null(map_to_label)){
    # 1) Keep only the mapped value that are valid
    map.cols <- c(map.AB.cols, map.A.cols, map.B.cols)
    data_map_tmp <- data_map
    data_map_tmp[data_map$map.type == "ND", c(map.AB.cols, map.A.cols, map.B.cols)] <- "ND"
    data_map_tmp[data_map$map.type == "strict", c(map.A.cols, map.B.cols)] <- "ND"
    data_map_tmp[data_map$map.type == "loose_A", c(map.AB.cols, map.B.cols)] <- "ND"
    data_map_tmp[data_map$map.type == "loose_B", c(map.AB.cols, map.A.cols)] <- "ND"

    # 2) Mapping
    factor.cols <- colnames(map_to_label)[sapply(map_to_label, is.factor)]
    other.cols <- colnames(map_to_label)[!sapply(map_to_label, is.factor)]
    for (col_ in factor.cols){
      data_map[[paste0(col_, "_all")]] <- apply(
        X = data_map_tmp[,map.cols], MARGIN = 1,
        FUN = .func_map_val_to_label_factor, map_to_label = map_to_label, col_label = col_)
      data_map[[col_]] <- gsub("[|].*", "", data_map[[paste0(col_, "_all")]])
      data_map[[col_]][data_map[[col_]]==""] <- as.character(NA)
    }
    for (col_ in other.cols){
      data_map[[col_]] <- apply(
        X = data_map_tmp[,map.cols], MARGIN = 1,
        FUN = .func_map_val_to_label, map_to_label = map_to_label, col_label = col_)
    }

    all.cols <- colnames(map_to_label)
    data_map[,paste0(all.cols, "_strict")] <- data_map[,all.cols]
    data_map[data_map$map.type != "strict", paste0(all.cols, "_strict")] <- as.character(NA)
    data_map[,paste0(all.cols, "_loose")] <- data_map[,all.cols]
    data_map[data_map$map.type == "ND", paste0(all.cols, "_loose")] <- as.character(NA)

  }

  return(data_map)
}
