# Generate synthetic VDJ contig data for testing

Creates a small, realistic VDJ contig annotation data frame mimicking
10x Genomics format. Designed for unit testing with controllable edge
cases including dropout, heterozygous chains, and invalid VDJ sequences.

## Usage

``` r
generate_test_VDJ(
  n_cells = 100,
  n_clones = 10,
  samples = c("SampleA", "SampleB"),
  dropout_rate = 0.2,
  heterozygous_rate = 0.3,
  invalid_rate = 0.05,
  seed = 42
)
```

## Arguments

- n_cells:

  Number of cells to generate. Default is 50.

- n_clones:

  Number of distinct clonotypes. Default is 10.

- samples:

  Character vector of sample names. Default is
  `c("SampleA", "SampleB")`.

- dropout_rate:

  Proportion of cells with missing chain data (0-1). Default is 0.2.

- heterozygous_rate:

  Proportion of cells with heterozygous (dual allele) chains (0-1).
  Default is 0.3.

- invalid_rate:

  Proportion of contigs with invalid VDJ (missing V or J gene) (0-1).
  Default is 0.05.

- seed:

  Random seed for reproducibility. Default is 42.

## Value

A data frame in 10x Genomics VDJ contig annotation format with columns:
barcode, is_cell, contig_id, high_confidence, length, chain, v_gene,
d_gene, j_gene, c_gene, full_length, productive, fwr1, fwr1_nt, cdr1,
cdr1_nt, fwr2, fwr2_nt, cdr2, cdr2_nt, fwr3, fwr3_nt, cdr3, cdr3_nt,
fwr4, fwr4_nt, reads, umis, raw_clonotype_id, raw_consensus_id,
exact_subclonotype_id, origin.

## Examples

``` r
test_data <- generate_test_VDJ(n_cells = 20, n_clones = 5, seed = 123)
head(test_data)
#>      barcode is_cell           contig_id high_confidence length chain v_gene
#> 1 CELL0001-1    true CELL0001-1_contig_1            true    514   TRB TRBV21
#> 2 CELL0002-1    true CELL0002-1_contig_2            true    604   TRA  TRAV3
#> 3 CELL0002-1    true CELL0002-1_contig_3            true    417   TRB  TRBV9
#> 4 CELL0003-1    true CELL0003-1_contig_4            true    507   TRA  TRAV3
#> 5 CELL0003-1    true CELL0003-1_contig_5            true    698   TRB  TRBV9
#> 6 CELL0004-1    true CELL0004-1_contig_6            true    564   TRA TRAV14
#>   d_gene j_gene c_gene full_length productive                        fwr1
#> 1  TRBD2  TRBJ1  TRBC2        true       true   CMAWVSDACCIYQPHIDSANVQIFF
#> 2   <NA>  TRAJ6   TRAC        true      false   CFPNIGSMNSRANLLEKYFKWDCNF
#> 3  TRBD1  TRBJ6  TRBC1       false       true      CVHHNKSDDQCEALGDNANPVF
#> 4   <NA>  TRAJ6   TRAC        true       true CNPGLILGWTCWWYPAGFVHTRFFPIF
#> 5  TRBD2  TRBJ6  TRBC2        true       true      CDQLYRHSNYKDWPITELTLFF
#> 6   <NA> TRAJ14   TRAC        true       true   CGDYQKYEHWIGCVHDDPIKNWCNF
#>                                                        fwr1_nt       cdr1
#> 1 CATGAAGGTCACGGGCAGGATATCTAGATTATACACTTCTCTCGACAACAGCTAAGGAAC   CSHEPNPF
#> 2 GAGAGAATTAAAGAATAAGTCCGTCGACTCTTTTCAGTGTAGTCGTGGTCTCCTGACAAC  CPVRPQVNF
#> 3 CACAGCTCGCGGTAGTAGTGGCGGCCTCGCGGTAACCTTGAAACAAGCTGTAACTAAGTC    CCLCGPF
#> 4 GACCTCCCCATGAGCTTCGCGCTGTCAGTTTCTTGCTGACCCCAAGGGCTTCCCCCCCGT    CTSFDGF
#> 5 GACACTCCTCAGCGGAGCATATACGTTCTGGCCTGTGAGAATAATCAGCTTGTAGGACCA  CPKTRQTTF
#> 6 GATTGTCGGGTCGTCGGCACACTTGCCTCCCTTTTTGGACAAGAGAAAATCTTTCAGGGT CEQITCFDFF
#>                    cdr1_nt                  fwr2
#> 1 CTTGACCGAGGAGGGGGTGTGCTA CDNQFQGKAYQWFCAYEMFEF
#> 2 TGTCCCAGAGGTCGGTCGATTAAG   CDRDCFFHCSIVWKWQWFF
#> 3 CGCGATGTGCACCGCTGGGAATCG CTGFTSKMFQFMPYEEMYGLF
#> 4 CAGGGCATCGGAGGCACAACTAGC     CEWKCTNQLHSCLWCIF
#> 5 GAACACCCACATCGGACCAGTCAG    CYKKTKLFMSQQPRQHLF
#> 6 CAGCTCACGGCTCTGCCGATAAGC     CWHISNLIYVYGHCTNF
#>                                              fwr2_nt       cdr2
#> 1 GGAGAAGAAGCAAACCCTTTCTGGTCACGAGGGACTGGCCAAACGGCCCA    CGEDMLF
#> 2 ACGCCCGTCGGGACAGAAGTTATCGAAACACAATTTCGTCTTATCGCAGG  CFIRQSAQF
#> 3 CGAGCAACTGGATCTGCCCTCCTAAAAACAGAGAAAGTAGGAAGGTTAAC  CLSAGRIAF
#> 4 AAGCTCGCCGTTTCCAGGCGGCTTACTGAGTGTAAATATAGAGCGCTTTC  CKREKFQTF
#> 5 ACAATCATATCGGGCAATATAAATTCCCGCGCGAATACAACTAAAAGTTA CGELRCHVKF
#> 6 ACTCTCGTGATCCTCAGTGGCACTTAATCGGGCTTGCAGAGGGCCTCACA CRKQIQTTHF
#>                    cdr2_nt                                      fwr3
#> 1 CAGACCAGTGGCGTAGCCGTTCCC CHHGTYAFLPLCRRIDCGQRWKRSIKSKCKMMEHVMAMGLF
#> 2 GGGAACGGAGCCTGCGAAAATTGG     CASIFLKFQWASGYIQQIMYFQAAFDEDQVLNSDNFF
#> 3 TATTATAGTGACAGTTTGAGATTT   CETYGHHMAWITQHDARQNTIACFPWYKSATVLKMACMF
#> 4 TGACCGAGGTAGGCATTCTCGCGC    CPSCPPGWMNMAIAINISLRLIICQVMCWSCKFGVTNF
#> 5 AGCTTCCGGCGTCATAGCGTACAC  CDCDAIDGIRTYWVDSMFNDGELYKRLVMVCFFYWRCLYF
#> 6 AAACCACAAGATCGGGCTCTAGTA  CDEKMWKLSNPPREILFVQGWWVGIHNENTLILRLMVYNF
#>                                                                                                fwr3_nt
#> 1 ATGGTTCCACTGTGACTGACTGTAGGCCGACGTTTGAATTCGTTAGAGGGGTTCCCAACCCGTATTCTCAGGTTTAGATGCAATCGATCAAGCGGGGTGT
#> 2 ACAAGTCTAGACGTTTGAAATCCACTTTCGATCCCGAGAGTCATTATCATGTGGTTTAGAGTTCTTGTGGCGAATTAAGGATAATCTCTATTACGTATGC
#> 3 AATTGTGTCCCAAATACGTCGCGTGCGGGTTGCATTTCAAAATTCTTGTTAGGCATGTAAGAATTTGGTGCAAGTGTACATGGTCAAGTATTGTTTGCTC
#> 4 GTTTCGATGACGGTGCGTCAACGAGGTGCAAAGATATCACCTTTCGGTCTTCGAGCATGTTATCCGCTGGCTAATGTAGAAGACTCGCGGGGTAGTAAGT
#> 5 TGATGAGGGAAGTCGAACCCAAACGGCCGCGTGGAACTAGTAGAAGCCGTCACGGTCGCCGGATGCGTGAACGTTAAACCGCAAGTGTGGGACCTGAGGG
#> 6 AATTAAGGATGTGACTGGCTGTGTCCAAGGCTTTCTTACAGCGTTCGAGTCCAGTTGCAGATAACCACCAATCCCGTTCCTAACGTCGGTAATTGATACC
#>                  cdr3                                                   cdr3_nt
#> 1 CTDVCFNPLGLLGSSDEDF GTGTTGCTAGAGTCAAGATGCCCGCTTTGCGGTCAGGGGGGCCAACAAATATTCCCC
#> 2        CCEPFWYQDISF                      ATCGTTGATAGCCACCCGAACTGCCCGGGGCTGGTG
#> 3  CDSCRDKHKECNKMHPFF    CCTATGGGCAGGGGTGCAGTCAGACGTACTCGCGTTCCAGTGGAGGTGCGCGGA
#> 4        CCEPFWYQDISF                      TCTTTACCGGGTGTAACGGCCGTAAAAACCTATTGA
#> 5  CDSCRDKHKECNKMHPFF    ACGAGCGCGAAGATGAACACGGGCCTCATCTGGTGACCGACTCCATGAAAAGAT
#> 6 CNPVAGRKRSYGMIHSTVF TGCAGCTCTCCCCCGACTATCTGTGCGCTGCAACTTCGGGTCTATTATTCGAAACTG
#>             fwr4                        fwr4_nt reads umis raw_clonotype_id
#> 1     CMDTYDQFVF CTCCCTTTAAGAAGACCTGGTTACGGACGA   911   40       clonotype3
#> 2  CGHESEKFSMSKF GCGCGAAAAAGCTGCGGGATTCGGTAAAGC  4212   25       clonotype4
#> 3     CPRNKWWGPF AGATGCTAGTTCAACGCTACCACGGCACTA  4272   24       clonotype4
#> 4 CRGRCTDVLNIIKF CCGAGTGGGACACACATTACACTAGTAATA   450   43       clonotype4
#> 5 CMFACEGFWNALEF TGCCTCTTGCCTTCGCGGCTGATTTCGACA  3279    9       clonotype4
#> 6 CVIYYEKMSLEAMF TCTGATTAGGCCAGCAGGGTCCTCAGTGGC   386   15       clonotype3
#>         raw_consensus_id exact_subclonotype_id  origin
#> 1 clonotype3_consensus_1                 sub_3 SampleB
#> 2 clonotype4_consensus_1                 sub_4 SampleB
#> 3 clonotype4_consensus_1                 sub_4 SampleB
#> 4 clonotype4_consensus_1                 sub_4 SampleA
#> 5 clonotype4_consensus_1                 sub_4 SampleA
#> 6 clonotype3_consensus_1                 sub_3 SampleA
```
