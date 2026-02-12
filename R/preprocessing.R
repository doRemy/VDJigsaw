# Preprocessing functions: pivot_VDJ and validate_TCR
# @Author Rémy Pétremand

# ------------------------------------------------------------------------------
# Pivot VDJ data
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
#' @param cdr3_col Which CDR3 column to use for building the composite TCR
#'   chain identifier. Either \code{"cdr3_nt"} (nucleotide sequence, default)
#'   or \code{"cdr3"} (amino acid sequence).
#' @param verbose Logical. If \code{TRUE}, print progress messages.
#'
#' @return A data frame in wide format with one row per cell barcode. Contains
#'   columns for sample, barcode, TCR chain identifiers (TRA_1, TRA_2, TRB_1,
#'   TRB_2), contig metadata, and individual chain component columns.
#'
#' @export
pivot_VDJ <- function(VDJ_data, sample_col = NULL, is_cell = FALSE, high_confidence = FALSE, productive = FALSE, full_length = FALSE, cdr3_col = "cdr3_nt", verbose = TRUE){

  if (verbose) message("[VDJigsaw] Converting VDJ data to wide format...")

  if (!cdr3_col %in% c("cdr3", "cdr3_nt")) {
    stop("cdr3_col must be either 'cdr3' (amino acid) or 'cdr3_nt' (nucleotide). Got: '", cdr3_col, "'")
  }
  if (verbose) message("[VDJigsaw]   Using ", cdr3_col, " for CDR3 chain identifier")

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
  if (any(duplicated(paste(VDJ_data_wide$Sample, VDJ_data_wide$barcode, sep = "_")))){
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
    FUN = function(x) paste(c("v_gene", cdr3_col, "j_gene"), x, sep = "_"))
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

# ------------------------------------------------------------------------------
# Validate TCR
# ------------------------------------------------------------------------------

#' Validate and correct TCR chain data
#'
#' Validates TCR chain columns (TRA_1, TRA_2, TRB_1, TRB_2) by:
#' - Checking V-CDR3-J format validity
#' - Reordering chains into the canonical V_CDR3_J format,
#' - Removing duplicate alleles (where allele 1 equals allele 2)
#' - Swapping alleles when allele 2 is defined but allele 1 is missing.
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
