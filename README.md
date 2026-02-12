# VDJigsaw — Puzzling Together T-Cell Clones from Single-Cell TCR-Seq Data

VDJigsaw is a toolkit for processing single-cell TCR-sequencing data into clonotype groups with varying levels of stringency and to map them to various types of representation.

## Introduction

Single-cell TCR sequencing (scTCR-seq) is a powerful approach for exploring T-cell receptor diversity and identifying clonally related T cells in humans. 
A common challenge is that a single cell barcode may be associated with one or multiple TCR alpha (TRA) and/or beta (TRB) chains, since some T cells are expected to carry heterozygous TCR-alpha and/or TCR-beta chains. 
Moreover, due to single-cell sequencing dropout and low TCR-chain detection sensitivity, chain information is frequently missing for individual cells. 
This makes clonotype assignment non-trivial.

For instance, consider a cell for which only one TRA chain allele and both alleles of the TRB chain were captured. 
Since no second TRA allele is present, two assumptions are possible:

- **Assumption 1 (homozygous):** The TCR-alpha chain is homozygous, and the captured allele is the only one.
- **Assumption 2 (heterozygous):** The TCR-alpha chain is heterozygous, and the second allele was simply not detected (dropout).

There is no way to distinguish between these two assumptions from a single cell in isolation. 
Current tools conservatively default to the first assumption and require an exact match on all chains or only consider the first alleles of both chains. 

However, we believe that this question can be addressed by leveraging the rest of the dataset. 
If we observe other cells that share the same three detected chains and additionally have a second TRA allele captured, it is reasonable to assume that all these cells belong to the same clone, and that the first cell simply failed to capture that chain (i.e., assumption 2). 
Otherwise, in the absence of such evidence, we default to assuming the second allele is identical to the first (i.e., assumption 1).

VDJigsaw addresses this by "puzzling" the different chains together into clone IDs at varying levels of stringency, allowing researchers to recover as many cells as possible into coherent clonotype definitions while controlling how permissive the matching is. 
We believe that this approach reduces the number of fragmented clonotypes originating from missing data and offers a more data-driven approach to clonotype assignment.

### Input Format

After pre-processing the scTCR-seq data, VDJigsaw processes the raw VDJ contigs data into a table where each row represents a cell (identified by its barcode) and the columns are:

| Column  | Description                        |
|---------|------------------------------------|
| `TRA_1` | Primary TCR alpha chain allele     |
| `TRA_2` | Secondary TCR alpha chain allele   |
| `TRB_1` | Primary TCR beta chain allele      |
| `TRB_2` | Secondary TCR beta chain allele    |

### Clonotype Stringency Levels

VDJigsaw groups cells into clonotypes using a hierarchy of stringency levels. Each successive level is more permissive, incorporating all the rules from the levels above it plus additional, more relaxed matching criteria. The goal is to treat missing values as potentially compatible rather than as mismatches, and to assign as many cells as possible to coherent clonotypes.

> **Note:** At any level of stringency, no clonotype will be defined by more than two TCR-alpha and two TCR-beta chains. This is enforced by grouping cells into a clonotype only when there are no conflicts (i.e., no two different values for the same chain position within a clonotype).

The stringency levels are named along two axes: **which chains are required** (both α and β, or just one) and **how many alleles per chain** (both alleles, or just one).

| Level | Name | Chains required | Alleles per chain |
|-------|------|----------------|-------------------|
| 1 | `dual_chain_dual_allele` | Both α and β | Both alleles (1 and 2) |
| 2 | `dual_chain_one_partial` | Both α and β | One chain can have a single allele |
| 3 | `dual_chain_both_partial` | Both α and β | Both chains can have a single allele |
| 4 | `single_chain_dual_allele` | α or β alone | Both alleles (1 and 2) |
| 5 | `single_chain_single_allele` | α or β alone | A single allele suffices |

#### `dual_chain_dual_allele` — Both chains, both alleles must match

The strictest level. All four chain positions must be identical between cells.

- `TRA_1`, `TRA_2`, `TRB_1`, and `TRB_2` must all match.

#### `dual_chain_one_partial` — Both chains required, one chain may have a missing allele

Relaxes the requirement on **one** of the two chains: either the alpha or the beta chain can have a missing secondary allele, as long as there are no conflicts.

- Full match on all four chains, **or**
- `TRA_1`, `TRB_1`, and `TRB_2` match with no conflicting `TRA_2` values (covers cases where `TRA_2` is missing), **or**
- `TRA_1`, `TRA_2`, and `TRB_1` match with no conflicting `TRB_2` values (covers cases where `TRB_2` is missing).

#### `dual_chain_both_partial` — Both chains required, both may have a missing allele

Relaxes the requirement on **both** chains: the secondary allele can be missing on both the alpha and beta chain simultaneously.

- All rules from `dual_chain_one_partial`, **or**
- `TRA_1` and `TRB_1` match with no conflicts on secondary chains.

#### `single_chain_dual_allele` — Only one chain required, but both alleles must match

No longer requires both receptor chains. A clonotype can be defined by a single chain type (alpha or beta), as long as both alleles of that chain match.

- All rules from `dual_chain_both_partial`, **or**
- `TRA_1` and `TRA_2` match with no conflicts on beta chains, **or**
- `TRB_1` and `TRB_2` match with no conflicts on alpha chains.

#### `single_chain_single_allele` — The most permissive level

A clonotype can be defined by a single allele of a single chain. This maximizes cell recovery but should be used with caution.

- All rules from `single_chain_dual_allele`, **or**
- `TRA_1` alone matches with no conflicts on other chains, **or**
- `TRB_1` alone matches with no conflicts on other chains.

## Installation

You can install VDJigsaw directly from GitHub: 

```r
# install.packages("devtools")
devtools::install_github("doRemy/VDJigsaw")
```

## Getting Started

```r
library(VDJigsaw)

# 1. Load VDJ filtered contigs
VDJ.contigs <- read.csv("filtered_contig_annotations.csv")

# 2. Assign clonotypes
VDJigsaw_res <- assign_clonotype(
  VDJ_data = VDJ.contigs,
  sample_col = "orig.ident",
  verbose = TRUE)

# 3. Retrieve outputs
TCR_data <- VDJigsaw_res$TCR_data
reference_tables <- VDJigsaw_res$ref_tables

# 4. Visualize
plot_stringency_summary(TCR_data)
plot_clone_composition(reference_tables$dual_chain_dual_allele)
plot_clone_composition(reference_tables$single_chain_single_allele)
```

For a full walkthrough with real data, see the [Getting Started vignette](vignettes/VDJigsaw_Vignette_01.Rmd).

## Output

`assign_clonotype()` returns a list with two elements:

- **`TCR_data`** — A data frame with one row per cell barcode, containing the TCR alpha/beta chains and a clone ID column for each stringency level.
- **`ref_tables`** — A named list of reference tables (one per stringency level), each mapping a `CloneID` to its TCR chain composition.

## Method

For each stringency level, VDJigsaw builds a pairwise similarity matrix between cells based on matching chain criteria (see [Stringency Levels](#clonotype-stringency-levels) above). This matrix is used to construct a graph where edges connect compatible cells. Connected components in this graph define the clonotypes.

Critically, any group containing a conflict — more than two distinct alleles for the same chain position — is not considered a valid match and is split accordingly. This ensures that no clonotype is defined by more than two TCR-alpha and two TCR-beta chains.

## Limitations

The main limitation of this approach is **doublets**. When two cells are captured in the same droplet, their TCR chains get mixed, potentially joining unrelated chains into a single clonotype. We recommend running VDJigsaw after quality control and doublet removal on your scRNA-seq data.

## Citation

If you use VDJigsaw in your research, please cite this repository:

> Pétremand, R. VDJigsaw: Puzzling Together T-Cell Clones from Single-Cell TCR-Seq Data. https://github.com/doRemy/VDJigsaw

## License

This project is licensed under the [GPL-3.0 License](LICENSE).

## Contributing

VDJigsaw is not currently accepting external contributions. This may change in the future — stay tuned.

## Transparency Note

AI (Claude, Anthropic) was used to assist with documentation writing and unit testing. All scientific content, methodology, and core implementation are the work of the author.
