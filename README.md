<div align="center">

# LaScaMol MR

[![Build Status](https://github.com/SamuelMathieu-code/MrPainter.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/SamuelMathieu-code/MrPainter.jl/actions/workflows/CI.yml?query=branch%3Amain)

</div>

## Overview

LaScaMolMR.jl (Large Scale Molecular Mendelian Randomization) is a threaded Mendelian Randomization (MR) package that is focused on the generation of transcriptome wide / molecular MR analysies. Although it provides interface for most common MR regression estimators (Inverse Variance Weighted, Weighted Median, Egger, Wald), its intended use is to enable fast Omic-wide Mendelian Randomization studies. The rise of large genetic cohort data has benefited the statistical power of Genome Wide Association Studies (GWAS) and Quantitative Trait Loci (QTL). Thus enabling findings in extensive studies such as Transcriptome Wide MR (TWMR), or mediation analyses between different levels of phenotypes. LaScaMolMR.jl provides a fast and efficient framework (still under developpement) to such analyses, allowing users to test different parameters.

## Example


For a QTL dataset composed as follows :

```
base/folder
├── exposureA.txt
├── exposureB.txt
└── exposureC.txt
```

With example data composed like this (tab separated) :

```
chr	pos	A1	A2	beta	se	some_useless_column pval
1	10511	A	G	-0.176656	0.136297	.   0.194939
1	10642	A	G	-0.724554	0.345390	..  0.035924
1	11008	G	C	-0.017786	0.016673	... 0.286088
1	11012	G	C	-0.017786	0.016673	..  0.286088
1	13110	A	G	0.013272	0.021949	.   0.545411
1	13116	G	T	-0.027802	0.013111	..  0.0339672
1	13118	G	A	-0.027802	0.013111	... 0.0339672
1	13259	A	G	-0.122207	0.210776	..  0.562052
1	13273	C	G	0.007077	0.015337	.   0.644463
```

and a GWAS of outcome composed similarly but comma separated, the following code will generate a trans-MR study with default parameters :

```julia
using LaScaMolMR
```

1. Describe exposure data.

```julia
path_pattern = ["exposure", TRAIT_NAME, ".txt"]
columns = Dict(1 => CHR, 2 => POS, 3 => A_EFFECT, 4 => A_OTHER, 5 => BETA, 6 => SE, 8 => PVAL)

trait_v = ["A", "B", "C"] # Chromosome and TSS informations are not relevant in Trans setting.

exposure::QTLStudy = QTLStudy_from_pattern("base/folder/", 
                                            path_pattern, 
                                            trait_v, chr_v = nothing, 
                                            tss_v = nothing, 
                                            columns, separator = '\t', 
                                            only_corresp_chr = false)
```

2. Describe outcome data.

```julia
outcome = GWAS("/some/file", columns, separator = ',', trait_name = "Some Painful Disease")
```

3. Perform Medelian Randomization study by providing input formats and reference genotype data files.

```julia
# Plink 1.9 files base names for each chromosome (You can also use a single file)
plink_files = ["folder/basename_file_chr$(i)" for i in 1:22]

# Perform MR for every exposure - outcme pairs with default parameters
out_table = mrStudyTrans(exposure, outcome, plink_files)

# with MiLoP approach and other parameters :
out_table2 = mrStudyTrans(exposure, outcome, plink_files, 
                        approach = "MiLoP", 
                        r2_tresh = 0.01, 
                        p_tresh = 5e-8,
                        filter_beta_ratio = 1)
# Default p_tresh_MiLoP value is the same as p_tresh
```

### MiLoP approach

Mitigated Local Pleiotropy (MiLoP) approach modifies the potential IV selection process to remove instrument variables associated to more than 1 exposure at `p_tresh_MiLoP` significance level.


## Enhancement ideas

- Multivariate MR & TWMR according to [Porcu et al.](https://pubmed.ncbi.nlm.nih.gov/31341166/)
- Mediation analysis inspired by [Auwerx et al.](https://elifesciences.org/articles/81097)
- Web interface for locus visualization with [Genie.jl](https://genieframework.com)
