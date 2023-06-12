<div align="center">

# LaScaMol MR

[![Build Status](https://github.com/SamuelMathieu-code/MrPainter.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/SamuelMathieu-code/MrPainter.jl/actions/workflows/CI.yml?query=branch%3Amain)

</div>

## Overview

LaScaMol.jl is a distrbuted Mendelian Randomization (MR) package that is focused on the generation of transcriptome wide / molecular MR analysies.

<span style="color:green">
coded 
</span> :

- inputs : modular inputs of GWAS and QTL studies.
- ld : computation of the ld composite
- MrPref : implementation of heterogeneity tests and iv-regression methods.

<span style="color:yellow">
ungoing 
</span> :

- MrStudyCis : Input parallel reading and filtering
- NaiveCis : MR study : clumping + MR

<span style="color:red">
planned 
</span> :

- StrictCis : Strict approach to MR Study to handle local internal pleiotropy
- SecondChanceCis : alternative aproach to handle local internal pleiotropy

<span style="color:purple">
enventualities
</span> :

- TransMR : trans selection of ivs
- ReverseMR : Reverse causality investigation
- PlotsMR : Graphical illustration of the results of an MR study