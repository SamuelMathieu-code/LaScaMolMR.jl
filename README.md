<div align="center">

# MrPainter

[![Build Status](https://github.com/SamuelMathieu-code/MrPainter.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/SamuelMathieu-code/MrPainter.jl/actions/workflows/CI.yml?query=branch%3Amain)

</div>

### Overview

MrPainter is a distrbuted Mendelian Randomization (MR) package that is focused on the generation of molecular MR analysies.

four parts :
    - inputs : modular inputs of GWAS and QTL studies.
    - ivSelect : selection process of ivs: cis, trans (& others in the future?)
    - MrPref : implementation of heterogeneity tests and iv-regression methods.
    - MrStudy : Large scale parallelisation of MRs from QTL to GWAS and GWAS to QTL.