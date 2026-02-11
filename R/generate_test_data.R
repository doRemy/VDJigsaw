#' Generate synthetic VDJ contig data for testing
#'
#' Creates a small, realistic VDJ contig annotation data frame mimicking 10x
#' Genomics format. Designed for unit testing with controllable edge cases
#' including dropout, heterozygous chains, and invalid VDJ sequences.
#'
#' @param n_cells Number of cells to generate. Default is 50.
#' @param n_clones Number of distinct clonotypes. Default is 10.
#' @param samples Character vector of sample names. Default is
#'   \code{c("SampleA", "SampleB")}.
#' @param dropout_rate Proportion of cells with missing chain data (0-1).
#'   Default is 0.2.
#' @param heterozygous_rate Proportion of cells with heterozygous (dual allele)
#'   chains (0-1). Default is 0.3.
#' @param invalid_rate Proportion of contigs with invalid VDJ (missing V or J
#'   gene) (0-1). Default is 0.05.
#' @param seed Random seed for reproducibility. Default is 42.
#'
#' @return A data frame in 10x Genomics VDJ contig annotation format with
#'   columns: barcode, is_cell, contig_id, high_confidence, length, chain,
#'   v_gene, d_gene, j_gene, c_gene, full_length, productive, fwr1, fwr1_nt,
#'   cdr1, cdr1_nt, fwr2, fwr2_nt, cdr2, cdr2_nt, fwr3, fwr3_nt, cdr3,
#'   cdr3_nt, fwr4, fwr4_nt, reads, umis, raw_clonotype_id, raw_consensus_id,
#'   exact_subclonotype_id, origin.
#'
#' @examples
#' test_data <- generate_test_VDJ(n_cells = 20, n_clones = 5, seed = 123)
#' head(test_data)
#'
#' @export
generate_test_VDJ <- function(n_cells = 50, n_clones = 10,
                               samples = c("SampleA", "SampleB"),
                               dropout_rate = 0.2,
                               heterozygous_rate = 0.3,
                               invalid_rate = 0.05,
                               seed = 42) {
  set.seed(seed)

  # Define realistic gene pools
  tra_v_genes <- paste0("TRAV", 1:20)
  tra_j_genes <- paste0("TRAJ", 1:15)
  trb_v_genes <- paste0("TRBV", 1:25)
  trb_d_genes <- paste0("TRBD", 1:2)
  trb_j_genes <- paste0("TRBJ", 1:10)
  tra_c_genes <- c("TRAC")
  trb_c_genes <- c("TRBC1", "TRBC2")

  # Generate random CDR3 amino acid sequences
  aa_pool <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
               "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  random_cdr3 <- function(n = 1, len_range = c(10, 18)) {
    vapply(seq_len(n), function(i) {
      len <- sample(len_range[1]:len_range[2], 1)
      paste0("C", paste(sample(aa_pool, len, replace = TRUE), collapse = ""), "F")
    }, character(1))
  }

  random_nt <- function(n = 1, len = 30) {
    vapply(seq_len(n), function(i) {
      paste(sample(c("A", "T", "G", "C"), len, replace = TRUE), collapse = "")
    }, character(1))
  }

  # Build clone definitions: each clone has a fixed set of chains
  clones <- data.frame(
    clone_id = seq_len(n_clones),
    tra1_v = sample(tra_v_genes, n_clones, replace = TRUE),
    tra1_cdr3 = random_cdr3(n_clones),
    tra1_j = sample(tra_j_genes, n_clones, replace = TRUE),
    tra2_v = sample(tra_v_genes, n_clones, replace = TRUE),
    tra2_cdr3 = random_cdr3(n_clones),
    tra2_j = sample(tra_j_genes, n_clones, replace = TRUE),
    trb1_v = sample(trb_v_genes, n_clones, replace = TRUE),
    trb1_cdr3 = random_cdr3(n_clones),
    trb1_j = sample(trb_j_genes, n_clones, replace = TRUE),
    trb2_v = sample(trb_v_genes, n_clones, replace = TRUE),
    trb2_cdr3 = random_cdr3(n_clones),
    trb2_j = sample(trb_j_genes, n_clones, replace = TRUE),
    stringsAsFactors = FALSE
  )

  # Some clones are homozygous (TRA_1 == TRA_2 or TRB_1 == TRB_2)
  homozygous_idx <- sample(n_clones, ceiling(n_clones * 0.3))
  clones$tra2_v[homozygous_idx] <- clones$tra1_v[homozygous_idx]
  clones$tra2_cdr3[homozygous_idx] <- clones$tra1_cdr3[homozygous_idx]
  clones$tra2_j[homozygous_idx] <- clones$tra1_j[homozygous_idx]

  # Assign cells to clones and samples
  cell_clone <- sample(n_clones, n_cells, replace = TRUE)
  cell_sample <- sample(samples, n_cells, replace = TRUE)

  # Decide cell types
  is_heterozygous_tra <- runif(n_cells) < heterozygous_rate
  is_heterozygous_trb <- runif(n_cells) < heterozygous_rate
  is_dropout <- runif(n_cells) < dropout_rate

  # Build contigs
  contigs <- list()
  contig_counter <- 0

  for (i in seq_len(n_cells)) {
    bc <- sprintf("CELL%04d-1", i)
    cl <- clones[cell_clone[i], ]
    smp <- cell_sample[i]
    clonotype_id <- paste0("clonotype", cell_clone[i])

    # Decide which chains to include
    include_tra1 <- TRUE
    include_tra2 <- is_heterozygous_tra[i]
    include_trb1 <- TRUE
    include_trb2 <- is_heterozygous_trb[i]

    # Apply dropout (remove some chains)
    if (is_dropout[i]) {
      dropout_type <- sample(c("no_tra", "no_trb", "no_tra2", "no_trb2"), 1)
      if (dropout_type == "no_tra") {
        include_tra1 <- FALSE
        include_tra2 <- FALSE
      } else if (dropout_type == "no_trb") {
        include_trb1 <- FALSE
        include_trb2 <- FALSE
      } else if (dropout_type == "no_tra2") {
        include_tra2 <- FALSE
      } else if (dropout_type == "no_trb2") {
        include_trb2 <- FALSE
      }
    }

    chains_to_add <- list()
    if (include_tra1) {
      chains_to_add[[length(chains_to_add) + 1]] <- list(
        chain = "TRA", v = cl$tra1_v, cdr3 = cl$tra1_cdr3, j = cl$tra1_j,
        d = NA, c_gene = "TRAC", umis = sample(5:50, 1)
      )
    }
    if (include_tra2) {
      chains_to_add[[length(chains_to_add) + 1]] <- list(
        chain = "TRA", v = cl$tra2_v, cdr3 = cl$tra2_cdr3, j = cl$tra2_j,
        d = NA, c_gene = "TRAC", umis = sample(1:20, 1)
      )
    }
    if (include_trb1) {
      chains_to_add[[length(chains_to_add) + 1]] <- list(
        chain = "TRB", v = cl$trb1_v, cdr3 = cl$trb1_cdr3, j = cl$trb1_j,
        d = sample(trb_d_genes, 1), c_gene = sample(trb_c_genes, 1), umis = sample(5:50, 1)
      )
    }
    if (include_trb2) {
      chains_to_add[[length(chains_to_add) + 1]] <- list(
        chain = "TRB", v = cl$trb2_v, cdr3 = cl$trb2_cdr3, j = cl$trb2_j,
        d = sample(trb_d_genes, 1), c_gene = sample(trb_c_genes, 1), umis = sample(1:20, 1)
      )
    }

    # Skip cells with no chains at all
    if (length(chains_to_add) == 0) next

    for (ch in chains_to_add) {
      contig_counter <- contig_counter + 1

      # Optionally make some contigs invalid (missing v or j gene)
      v_gene <- ch$v
      j_gene <- ch$j
      if (runif(1) < invalid_rate) {
        if (runif(1) < 0.5) v_gene <- NA else j_gene <- NA
      }

      contigs[[contig_counter]] <- data.frame(
        barcode = bc,
        is_cell = "true",
        contig_id = sprintf("%s_contig_%d", bc, contig_counter),
        high_confidence = sample(c("true", "false"), 1, prob = c(0.95, 0.05)),
        length = sample(400:700, 1),
        chain = ch$chain,
        v_gene = v_gene,
        d_gene = ch$d,
        j_gene = j_gene,
        c_gene = ch$c_gene,
        full_length = sample(c("true", "false"), 1, prob = c(0.9, 0.1)),
        productive = sample(c("true", "false"), 1, prob = c(0.95, 0.05)),
        fwr1 = random_cdr3(1, c(20, 25)),
        fwr1_nt = random_nt(1, 60),
        cdr1 = random_cdr3(1, c(5, 8)),
        cdr1_nt = random_nt(1, 24),
        fwr2 = random_cdr3(1, c(15, 20)),
        fwr2_nt = random_nt(1, 50),
        cdr2 = random_cdr3(1, c(5, 8)),
        cdr2_nt = random_nt(1, 24),
        fwr3 = random_cdr3(1, c(30, 40)),
        fwr3_nt = random_nt(1, 100),
        cdr3 = ch$cdr3,
        cdr3_nt = random_nt(1, nchar(ch$cdr3) * 3),
        fwr4 = random_cdr3(1, c(8, 12)),
        fwr4_nt = random_nt(1, 30),
        reads = sample(100:5000, 1),
        umis = ch$umis,
        raw_clonotype_id = clonotype_id,
        raw_consensus_id = paste0(clonotype_id, "_consensus_1"),
        exact_subclonotype_id = paste0("sub_", cell_clone[i]),
        origin = smp,
        stringsAsFactors = FALSE
      )
    }
  }

  result <- do.call(rbind, contigs)
  rownames(result) <- NULL
  return(result)
}
