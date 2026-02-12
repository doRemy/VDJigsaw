# Validate and correct TCR chain data

Validates TCR chain columns (TRA_1, TRA_2, TRB_1, TRB_2) by:

- Checking V-CDR3-J format validity

- Reordering chains into the canonical V_CDR3_J format,

- Removing duplicate alleles (where allele 1 equals allele 2)

- Swapping alleles when allele 2 is defined but allele 1 is missing.

## Usage

``` r
validate_TCR(TCR_data, remove_invalid_VDJ = FALSE, verbose = TRUE)
```

## Arguments

- TCR_data:

  A data frame with columns TRA_1, TRA_2, TRB_1, TRB_2 containing TCR
  chain identifiers in V_CDR3_J format.

- remove_invalid_VDJ:

  Logical. If `TRUE`, set invalid VDJ sequences (those with missing V,
  CDR3, or J components) to `NA`. If `FALSE`, keep them as-is and flag
  them.

- verbose:

  Logical. If `TRUE`, print progress messages.

## Value

The input data frame with corrected TCR chain columns and additional
logical columns (TRA_1.invalid, TRA_2.invalid, TRB_1.invalid,
TRB_2.invalid) indicating which chains had invalid VDJ format.
