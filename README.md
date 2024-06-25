<div align="center">

# LaScaMol MR

[![Build Status](https://github.com/SamuelMathieu-code/MrPainter.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/SamuelMathieu-code/MrPainter.jl/actions/workflows/CI.yml?query=branch%3Amain)

</div>

## Overview

LaScaMolMR.jl (Large Scale Molecular Mendelian Randomization) is a threaded Mendelian Randomization (MR) package that is focused on the generation of transcriptome wide / molecular MR analysies. Although it provides interface for most common MR regression estimators (Inverse Variance Weighted, Weighted Median, Egger, Wald), its intended use is to enable fast Omic-wide Mendelian Randomization studies. The rise of large genetic cohort data has benefited the statistical power of Genome Wide Association Studies (GWAS) and Quantitative Trait Loci (QTL). Thus enabling findings in extensive studies such as Transcriptome Wide MR (TWMR), or mediation analyses between different levels of phenotypes ([Porcu et al.](https://elifesciences.org/articles/81097)). LaScaMolMR.jl provides a fast and efficient framework (still under developpement) to such analyses, allowing users to test different parameters.

## Tutorial

### Example 1 : Cis-MR

For a QTL dataset composed as follows :

```
base/folder
├── exposureA
│   ├── exposureA_chr1.txt
│   └── exposureA_chr2.txt
├── exposureB
│   ├── exposureB_chr1.txt
│   └── exposureB_chr2.txt
└── exposureC
    ├── exposureC_chr3.txt
    └── exposureC_chr4.txt
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

and a GWAS of outcome composed similarly but comma separated, the following code will generate a cis-MR study with default parameters :

```julia
using LaScaMolMR

# Describing the exposure and the outcome:

# Description of pattern of information composing the path to each file.
path_pattern = ["exposure", TRAIT_NAME, "/exposure", TRAIT_NAME, "_chr", CHR, ".txt"]

# Description of information contained in columns of QTL files. 
# See GenVarInfo in documentation for more details.
columns = Dict(1 => CHR, 2 => POS, 3 => A_EFFECT, 4 => A_OTHER, 5 => BETA, 6 => SE, 8 => PVAL)

# compose of load traits to analyse and their TSS. Using a subset of possible exposure will 
# use only corresponding files. We suggest DataFrames.jl + CSV.jl or DelimitedFiles.jl 
# packages to load these from a file.
trait_v = ["A", "B", "C"] # exposure trait identifiers (vector)
chr_v = [1, 2, 3] # chromosome corresponding to exposure protein or gene (vector) 
                  # (only the corresponding file will be used)
tss_v = [45287, 984276, 485327765] # position of TSS (vector)

# QTLStudy type for exposure and GWAS for outcome
exposure::QTLStudy = QTLStudy_from_pattern("base/folder/", 
                                            path_pattern, 
                                            trait_v, chr_v, tss_v, 
                                            columns, separator = '\t')

outcome = GWAS("/some/file", columns, separator = ',', trait_name = "Some Painful Desease")

# Plink 1.9 files base names for each chromosome (You can also use a single file)
plink_files = ["folder/basename_file_chr$(i)" for i in 1:22]

# Perform MR for every exposure - outcme pairs with default parameters
out_table = mrStudyCis(exposure, outcome, plink_files)

# with strict approach to internal pleiotropy and other parameters :
out_table2 = mrStudyCis(exposure, outcome, plink_files, 
                        approach = "strict", 
                        r2_tresh = 0.01, 
                        p_tresh = 5e-8, 
                        pval_bigfolat = true)
# This method provides many options related to the parameters, 
# rendered outputs and the MR methods used (see documentation for detailed information)

# You can also separate te IV filtration part of the Study from the clumping+MR part 
# by using the test and test-strict approaches :

potential_ivs = out_table = mrStudyCis(exposure, outcome, plink_files, approach = "test")

using SnpArrays

#load and format snp data from Plink 1.9
genotypes = [SnpData(SnpArrays.datadir(file)) for file in plink_files]
for snpdata in genotypes
    formatSnpData!(snpdata)
end

# Run clumping + MR analysis with pre-selected ivs and specified MR functions
table_mr_results = NaiveCis(potential_ivs, genotypes, mr_methods = [mr_ivw, mr_wald])
```

### Example 2 : Trans-MR

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

# Describing the exposure and the outcome:

# Description of pattern of information composing the path to each file.
path_pattern = ["exposure", TRAIT_NAME, ".txt"]

# Description of information contained in columns of QTL files. 
# See GenVarInfo in documentation for more details.
columns = Dict(1 => CHR, 2 => POS, 3 => A_EFFECT, 4 => A_OTHER, 5 => BETA, 6 => SE, 8 => PVAL)

# compose of load traits to analyse and their TSS. Using a subset of possible exposure 
# will use only corresponding files. We suggest DataFrames.jl + CSV.jl or DelimitedFiles.jl 
# packages to load these from a file.
trait_v = ["A", "B", "C"] # exposure trait identifiers (vector)

# QTLStudy type for exposure and GWAS for outcome

exposure::QTLStudy = QTLStudy_from_pattern("base/folder/", 
                                            path_pattern, 
                                            trait_v, chr_v = nothing, 
                                            tss_v = nothing, 
                                            columns, separator = '\t', 
                                            only_corresp_chr = false)

outcome = GWAS("/some/file", columns, separator = ',', trait_name = "Some Painful Desease")

# Plink 1.9 files base names for each chromosome (You can also use a single file)
plink_files = ["folder/basename_file_chr$(i)" for i in 1:22]

# Perform MR for every exposure - outcme pairs with default parameters
out_table = mrStudyTrans(exposure, outcome, plink_files)

# with strict approach to internal pleiotropy and other parameters :
out_table2 = mrStudyTrans(exposure, outcome, plink_files, 
                        approach = "strict", 
                        r2_tresh = 0.01, 
                        p_tresh = 5e-8, 
                        filter_beta_ratio = 1)
# This method provides many options related to the parameters, 
# rendered outputs and the MR methods used (see documentation for detailed information)
```

### `mr_*` regression functions

LaScaMolMR provides four MR regression functions : `mr_wald`, `mr_ivw`, `mr_egger` and `mr_wm` performing respectively the wald ratio of the first provided instrument, Inverse Variance Weighted, Egger and Weighted Median regressions.

Example :

```julia
out = mr_wm(beta_outcome, se_outcome, beta_exposure, se_exposure;
            iteration = 10000, seed = 42)
```

### mr_output standard functions :

The mr_output struct provies a standard interface to functions performing mendelian randomization. This allows for user to use its own MR functions for mrStudies. Any function receiving 4 vectors and having at least $\alpha$ as an option. Such a function should return an mr_output object.

Here is the definitioin of mr_output :

```julia
struct mr_output
    nivs::Int
    effect::Float64
    se_effect::Float64
    ci_low::Float64
    ci_high::Float64
    p::Float64
    intercept::Float64
    p_intercept::Float64
    ci_low_intercept::Float64
    ci_high_intercept::Float64
    heter_stat::Float64
    heter_p::Float64
end
```

Any function following this format will could be provided to NaiveCis/NaiveTrans/mrStudyCis/mrStudyTrans ;

```julia
function mr_something(beta_outcome, se_outcome, beta_exposure, se_exposure; 
                      α = 0.05, other_params = default_value
                     )::mr_output

    # Calculate values
    ...
    # end

    # Not calculated values can be replaced by NaN
    return mr_ouput(nivs, effect, ci_low, ci_high, p, NaN, NaN, NaN, NaN, heter_stat, heter_p)

end

output = NaiveCis(data, genotypes; mr_methods = [mr_ivw, mr_something])
```


## TODO

- Transformer strict en Milop -> double seuil.
- Simplifier et renommer NaiveCis et NaiveTrans ->  une seule fonction
- implémenter mvmr (seule)

