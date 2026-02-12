# VDJigsaw constants and package-level definitions
# @Author Rémy Pétremand

#' @importFrom dplyr %>% group_by group_by_at ungroup summarise mutate filter reframe if_else dense_rank n select all_of across where .data mutate_at n_distinct arrange desc count slice pull bind_rows transmute
#' @importFrom tidyr pivot_wider pivot_longer separate_wider_delim
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

# Get the clone-definition dataframe
.clone.definition.df <- rbind(
  c("ColName",              "CloneDef",                 "CloneDefRef",              "TRA_1", "TRA_2", "TRB_1", "TRB_2"),
  # dual_chain_dual_allele combinations
  c("Sample_A1_A2_B1_B2", "dual_chain_dual_allele",     "None",                     "TRA_1", "TRA_2", "TRB_1", "TRB_2"),
  c("Sample_A1_A2_B2_B1", "dual_chain_dual_allele",     "None",                     "TRA_1", "TRA_2", "TRB_2", "TRB_1"),
  c("Sample_A2_A1_B1_B2", "dual_chain_dual_allele",     "None",                     "TRA_2", "TRA_1", "TRB_1", "TRB_2"),
  c("Sample_A2_A1_B2_B1", "dual_chain_dual_allele",     "None",                     "TRA_2", "TRA_1", "TRB_2", "TRB_1"),
  # dual_chain_one_partial combinations built from dual_chain_dual_allele
  c("Sample_A1_A2_B1_NA", "dual_chain_one_partial",     "dual_chain_dual_allele",   "TRA_1", "TRA_2", "TRB_1", "NAcol"),
  c("Sample_A1_A2_B2_NA", "dual_chain_one_partial",     "dual_chain_dual_allele",   "TRA_1", "TRA_2", "TRB_2", "NAcol"),
  c("Sample_A2_A1_B1_NA", "dual_chain_one_partial",     "dual_chain_dual_allele",   "TRA_2", "TRA_1", "TRB_1", "NAcol"),
  c("Sample_A2_A1_B2_NA", "dual_chain_one_partial",     "dual_chain_dual_allele",   "TRA_2", "TRA_1", "TRB_2", "NAcol"),
  c("Sample_A1_NA_B1_B2", "dual_chain_one_partial",     "dual_chain_dual_allele",   "TRA_1", "NAcol", "TRB_1", "TRB_2"),
  c("Sample_A1_NA_B2_B1", "dual_chain_one_partial",     "dual_chain_dual_allele",   "TRA_1", "NAcol", "TRB_2", "TRB_1"),
  c("Sample_A2_NA_B1_B2", "dual_chain_one_partial",     "dual_chain_dual_allele",   "TRA_2", "NAcol", "TRB_1", "TRB_2"),
  c("Sample_A2_NA_B2_B1", "dual_chain_one_partial",     "dual_chain_dual_allele",   "TRA_2", "NAcol", "TRB_2", "TRB_1"),
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
.clone.definition.df <- as.data.frame(.clone.definition.df, stringsAsFactors = FALSE)
colnames(.clone.definition.df) <- .clone.definition.df[1, ]
.clone.definition.df <- .clone.definition.df[-1, ]

# Get the clone mapping def order
.mapping.levels <- c(
  "dual_chain_dual_allele", 
  "dual_chain_one_partial", 
  "dual_chain_both_partial", 
  "single_chain_dual_allele", 
  "single_chain_single_allele")

# Set CloneDef columns as factor:
.clone.definition.df$CloneDef <- factor(.clone.definition.df$CloneDef, levels = .mapping.levels)
