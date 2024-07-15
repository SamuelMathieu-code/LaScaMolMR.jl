# LaScaMolMR

## Overview 

LaScaMolMR.jl (Large Scale Molecular Mendelian Randomization) is a threaded Mendelian Randomization (MR) package that is focused on the generation of transcriptome wide / molecular MR analysies. Although it provides an interface for most common MR regression estimators (Inverse Variance Weighted, Weighted Median, Egger, Wald), its intended use is to enable fast Omic-wide Mendelian Randomization studies. The rise of large genetic cohort data has benefited the statistical power of Genome Wide Association Studies (GWAS) and Quantitative Trait Loci (QTL). Thus enabling findings in extensive studies such as Transcriptome Wide MR (TWMR). LaScaMolMR.jl provides a fast and efficient framework to such analyses, allowing users to customize the parameters of the study.

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
chr	pos	effect_allele	other_allele	beta	se	some_other_column pval
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
```

1. Describe the exposure data structure (files and columns).

```julia
path_pattern = ["exposure", TRAIT_NAME, "/exposure", TRAIT_NAME, "_chr", CHR, ".txt"]
columns = Dict(1 => CHR, 2 => POS, 3 => A_EFFECT, 4 => A_OTHER, 5 => BETA, 6 => SE, 8 => PVAL)

trait_v = ["A", "B", "C"] # exposure trait identifiers
chr_v = [1, 2, 3] # chromosome corresponding to exposure protein or gene
tss_v = [45287, 984276, 485327765] # position of TSS, if relevant

exposure::QTLStudy = QTLStudy_from_pattern("base/folder/", 
                                            path_pattern, 
                                            trait_v, chr_v, tss_v, 
                                            columns, separator = '\t')
```

2. Describe the outcome data.

```julia
outcome = GWAS("/some/file", columns, separator = ',', trait_name = "Some Painful Disease")
```

3. Perform Medelian Randomization study by providing input formats and reference genotype data files.

```julia
# Plink 1.9 files base names for each chromosome (You can also use a single file)
plink_files = ["folder/basename_file_chr$(i)" for i in 1:22]

# Perform MR for every exposure - outcome pairs with default parameters
out_table = mrStudy(exposure, outcome, "cis", plink_files)

# with MiLoP approach to internal pleiotropy and other parameters :
out_table_MiLoP = mrStudy(exposure, outcome, "cis", plink_files, 
                        approach = "MiLoP", 
                        r2_tresh = 0.01, 
                        p_tresh = 5e-8,
                        p_tresh_MiLoP = 5e-4,
                        pval_bigfolat = true)
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
out_table = mrStudy(exposure, outcome, "trans", plink_files)

# with MiLoP approach and other parameters :
out_table2 = mrStudy(exposure, outcome, "trans", plink_files, 
                        approach = "MiLoP", 
                        r2_tresh = 0.01, 
                        p_tresh = 5e-8,
                        filter_beta_ratio = 1)
# Default p_tresh_MiLoP value is the same as p_tresh
```

### MiLoP approach

Mitigated Local Pleiotropy (MiLoP) approach modifies the potential IV selection process to remove instrument variables associated to more than 1 exposure at `p_tresh_MiLoP` significance level.

### `mr_*` regression functions

LaScaMolMR provides four MR regression functions : `mr_wald`, `mr_ivw`, `mr_egger` and `mr_wm` performing respectively the wald ratio of the first provided instrument, Inverse Variance Weighted, Egger and Weighted Median regressions.

Example :

```julia
out = mr_wm(beta_outcome, se_outcome, beta_exposure, se_exposure;
            iteration = 10000, seed = 42, α = 0.05)
```

### `mr_output` as a standard output for `mr_*` :

The mr_output struct provies a standard interface to functions performing mendelian randomization. This allows for user to use its own MR functions for `mrStudy`. Any function receiving 4 vectors, having at least $\alpha$ as an option and outputing a `mr_output` is valid.

Here is the definitioin of `mr_output` :

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

Any function following this format could be provided to [`clumpAndMR`](@ref)/[`mrStudy`](@ref) in the `mr_methods` option,
including user-defined functions.

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

## API

### Index

```@index
```

### Provide inputs

These functions allow users to specify the input formats. The QTL summary statistics can be united in a single file or divided by exposure, chromosome or both.

For the sake of simplicity, we consider a genetic association study as a GWAS when a single phenotype was studied and as a QTLStudy otherwise.

```@docs
GenVarInfo

GWAS

QTLStudy

QTLStudy_from_pattern
```

### Perform a MR Study

The function `mrStudy` allows the users to perform meta-analysis over many pairs of traits given the input format. The function has three implementations depending on the type of studies used as exposure and outcomes. (GWAS -> GWAS, QTL -> GWAS, GWAS -> QTL)

```@docs
mrStudy

mrStudyNFolds
```

### Linkage Desiquilibrium and MR

The `clumpAndMR` function performs clumping over every pair of (exposure, outcome) and calls Mendelian randomization fonctions. You can thus use it if potential IV selection was already performed.

```@docs
clumpAndMR
```

### Mendelian Randomization functions

```@docs
mr_output

mr_wald

mr_ivw

mr_egger

mr_wm
```

### Exported Utilities

#### Inputs

```@docs
nfolds
```
#### LD

```@docs
formatSnpData!

clump
```

## Contents

```@contents
```

### Authors 

**Samuel Mathieu**, Hippolyte Minvielle Moncla, Mewen Briend, Valentine Duclos, Anne Rufiange, Patrick Mathieu
