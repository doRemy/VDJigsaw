# VDJigsaw — Puzzling Together T-Cell Clones from Single-Cell TCR-Seq Data

VDJigsaw is a toolkit for processing single-cell TCR-sequencing data
into clonotype groups with varying levels of stringency and to map them
to various types of representation.

## Introduction

Single-cell T Cell Receptor sequencing (scTCR-seq) is a powerful
approach for exploring TCR diversity and identifying clonally related T
cells. TCRs are heterodimers composed of two protein chains: alpha and
beta. In most cases, a T cell is expected to express only one alpha-beta
pair. However, while it is unlikely that a cell expresses two distinct
beta chains, it is relatively common (~30%) for a cell to express two
alpha chains.

Clonotype assignment — grouping cells into the same clone — is a
challenging problem in scTCR-seq data. What makes it particularly
non-trivial is that chain information is frequently missing for
individual cells due to sequencing dropout and low TCR-chain detection
sensitivity.

As a concrete example: imagine two cells, Cell-1 and Cell-2, sharing the
same TCR-beta chain. Cell-1 has one TCR-alpha chain captured (alpha-1),
while Cell-2 has two (alpha-1 and alpha-2). Note that alpha-1 matches
between the two cells.

Since Cell-1 and Cell-2 share the same beta chain and at least one alpha
chain allele, the key question is whether:

- Cell-1 actually has a second alpha allele identical to alpha-1
  (homozygous), or
- Cell-1 is heterozygous with alpha-2 as its second allele, but alpha-2
  was not captured due to dropout.

Current tools conservatively default to the first assumption and require
an exact match on all chains, or only consider the primary alleles of
both chains.

However, given the immensely large diversity of TCRs generated through
VDJ recombination, it is reasonable to assume that Cell-1 and Cell-2
belong to the same clonotype — it is very unlikely that two independent
T-cell clones share the same beta chain and one of the same alpha
alleles by chance.

To test this, one could examine the other TCR sequences in the sample
for evidence refuting this hypothesis. If no conflicting observation
exists, it is reasonable to assume that these cells belong to the same
clone and that the missing chain in Cell-1 simply failed to be captured.

VDJigsaw addresses this by “puzzling” the different chains together into
clone IDs at varying levels of stringency, allowing researchers to
recover as many cells as possible into coherent clonotype definitions
while controlling how permissive the matching is.

This approach reduces the number of fragmented clonotypes caused by
missing data and offers a more data-driven method for clonotype
assignment.

See the [Method](#method) section below for more details.

## Installation

You can install VDJigsaw directly from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("doRemy/VDJigsaw")
```

## Getting Started

``` r
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

For a full walkthrough with real data, see the [Getting Started
vignette](https://doRemy.github.io/VDJigsaw/articles/VDJigsaw_Vignette_01.html).

## Output

[`assign_clonotype()`](https://doRemy.github.io/VDJigsaw/reference/assign_clonotype.md)
returns a list with two elements:

- **`TCR_data`** — A data frame with one row per cell barcode,
  containing the TCR alpha/beta chains and a clone ID column for each
  stringency level.
- **`ref_tables`** — A named list of reference tables (one per
  stringency level), each mapping a `CloneID` to its TCR chain
  composition.

## Method

### Input Format

After pre-processing the scTCR-seq data, VDJigsaw processes the raw VDJ
contigs data into a table where each row represents a cell (identified
by its barcode) and the columns are:

| Column  | Description                      |
|---------|----------------------------------|
| `TRA_1` | Primary TCR alpha chain allele   |
| `TRA_2` | Secondary TCR alpha chain allele |
| `TRB_1` | Primary TCR beta chain allele    |
| `TRB_2` | Secondary TCR beta chain allele  |

### Clonotype Stringency Levels

VDJigsaw groups cells into clonotypes using a hierarchy of stringency
levels. Each successive level is more permissive, incorporating all the
rules from the levels above it plus additional, more relaxed matching
criteria. The goal is to treat missing values as potentially compatible
rather than as mismatches, and to assign as many cells as possible to
coherent clonotypes.

> **Note:** At any level of stringency, no clonotype will be defined by
> more than two TCR-alpha and two TCR-beta chains. This is enforced by
> grouping cells into a clonotype only when there are no conflicts
> (i.e., no two different values for the same chain position within a
> clonotype).

The stringency levels are named along two axes: **which chains are
required** (both α and β, or just one) and **how many alleles per
chain** (both alleles, or just one).

| Level | Name                         | Chains required | Alleles per chain                    |
|-------|------------------------------|-----------------|--------------------------------------|
| 1     | `dual_chain_dual_allele`     | Both α and β    | Both alleles (1 and 2)               |
| 2     | `dual_chain_one_partial`     | Both α and β    | One chain can have a single allele   |
| 3     | `dual_chain_both_partial`    | Both α and β    | Both chains can have a single allele |
| 4     | `single_chain_dual_allele`   | α or β alone    | Both alleles (1 and 2)               |
| 5     | `single_chain_single_allele` | α or β alone    | A single allele suffices             |

#### `dual_chain_dual_allele` — Both chains, both alleles must match

The strictest level. All four chain positions must be identical between
cells.

- `TRA_1`, `TRA_2`, `TRB_1`, and `TRB_2` must all match.

#### `dual_chain_one_partial` — Both chains required, one chain may have a missing allele

Relaxes the requirement on **one** of the two chains: either the alpha
or the beta chain can have a missing secondary allele, as long as there
are no conflicts.

- Full match on all four chains, **or**
- `TRA_1`, `TRB_1`, and `TRB_2` match with no conflicting `TRA_2` values
  (covers cases where `TRA_2` is missing), **or**
- `TRA_1`, `TRA_2`, and `TRB_1` match with no conflicting `TRB_2` values
  (covers cases where `TRB_2` is missing).

#### `dual_chain_both_partial` — Both chains required, both may have a missing allele

Relaxes the requirement on **both** chains: the secondary allele can be
missing on both the alpha and beta chain simultaneously.

- All rules from `dual_chain_one_partial`, **or**
- `TRA_1` and `TRB_1` match with no conflicts on secondary chains.

#### `single_chain_dual_allele` — Only one chain required, but both alleles must match

No longer requires both receptor chains. A clonotype can be defined by a
single chain type (alpha or beta), as long as both alleles of that chain
match.

- All rules from `dual_chain_both_partial`, **or**
- `TRA_1` and `TRA_2` match with no conflicts on beta chains, **or**
- `TRB_1` and `TRB_2` match with no conflicts on alpha chains.

#### `single_chain_single_allele` — The most permissive level

A clonotype can be defined by a single allele of a single chain. This
maximizes cell recovery but should be used with caution.

- All rules from `single_chain_dual_allele`, **or**
- `TRA_1` alone matches with no conflicts on other chains, **or**
- `TRB_1` alone matches with no conflicts on other chains.

### Graph-Based Clonotype Assignment Algorithm

Clonotype assignment is performed **independently per sample** and
proceeds through the stringency levels **from strictest to most
permissive**. At each level, VDJigsaw:

1.  **Generates clone version strings** — For each cell, the relevant
    chain values (as defined by the current stringency level) are
    concatenated into a composite identifier. For example, at the
    strictest level, this is `Sample_TRA1_TRA2_TRB1_TRB2`; at more
    permissive levels, some chain positions are omitted.

2.  **Builds a pairwise similarity matrix** — Each cell’s full clone
    version string is compared against all other cells’ clone version
    strings for the current level. Two cells are considered compatible
    if any of their clone version strings match. Matches from stricter
    levels are carried forward, so compatibility is cumulative.

3.  **Constructs a graph** — The similarity matrix is used as an
    adjacency matrix to build an undirected graph (via
    [igraph](https://r.igraph.org/)), where nodes are cells and edges
    connect compatible cells.

4.  **Identifies connected components** — Each connected component in
    the graph defines a candidate clonotype. Cells with no edges (no
    match at any level) remain unassigned.

5.  **Resolves conflicts** — A connected component is checked for
    conflicts: if it contains more than two distinct alleles for any
    chain position (e.g., three different TRA_1 values), the component
    is invalid. In this case, VDJigsaw falls back to the clustering from
    the previous (stricter) stringency level for the affected cells,
    ensuring that no clonotype is defined by more than two alpha and two
    beta chain alleles.

This hierarchical approach means that cells first assigned at a strict
level retain that assignment, and each successive level only adds new
assignments for cells that were previously unmatched — while always
enforcing the two-allele-per-chain constraint.

## Limitations

The main limitation of this approach is **doublets**. When two cells are
captured in the same droplet, their TCR chains get mixed, potentially
joining unrelated chains into a single clonotype. We recommend running
VDJigsaw after quality control and doublet removal on your scRNA-seq
data.

## Citation

If you use VDJigsaw in your research, please cite this repository:

> Pétremand, R. VDJigsaw: Puzzling Together T-Cell Clones from
> Single-Cell TCR-Seq Data. <https://github.com/doRemy/VDJigsaw>

## License

This project is licensed under the [GPL-3.0
License](https://doRemy.github.io/VDJigsaw/LICENSE).

## Contributing

VDJigsaw is not currently accepting external contributions. This may
change in the future — stay tuned.

## Transparency Note

AI (Claude, Anthropic) was used to assist with documentation writing and
unit testing. All scientific content, methodology, and core
implementation are the work of the author.
