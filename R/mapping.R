# Post-pipeline mapping functions
# @Author Rémy Pétremand

# ------------------------------------------------------------------------------
# assign_clonotype_from_reference
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
  TCR_data  <- cbind(TCR_data,  .get_clone_versions(TCR_data = TCR_data))
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
  
  # Part 5: Return the mapped CloneID
  return(TCR_data$CloneID.map)
}

# ------------------------------------------------------------------------------
# annotate_by_clone
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# map_clonotypes_to_paired_TCR
# ------------------------------------------------------------------------------

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

  # Define some internal functions:
  .func_map_val <- function(x){
    x = sort(unique(x[x != "ND"]))
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

  # Note: order or operations matters. TCR-beta chain has priority in loose matching
  data_map$map.loose[mask.fill.mapA] <- data_map$map.A[mask.fill.mapA]
  data_map$map.loose[mask.fill.mapB] <- data_map$map.B[mask.fill.mapB]

  # Part 4: Get summary
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

  # Part 6: return the mapped data
  return(data_map)
}
