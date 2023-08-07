<div align="center">

# LaScaMol MR

[![Build Status](https://github.com/SamuelMathieu-code/MrPainter.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/SamuelMathieu-code/MrPainter.jl/actions/workflows/CI.yml?query=branch%3Amain)

</div>

## Overview

LaScaMolMR.jl is a threaded Mendelian Randomization (MR) package that is focused on the generation of transcriptome wide / molecular MR analysies.

<span style="color:green">
<b>coded </b>
</span> :

- inputs : modular inputs of GWAS and QTL studies.
- ld : computation of the ld composite
- MrPref : implementation of heterogeneity tests and iv-regression methods.
- mrStudyCis, mrStudyTrans, NaiveCis, NaiveTrans : Implemenation of IV selection in cis- and trans-. Supports both strict and naive aproaches for internal pleiotropy.

<span style="color:purple">
<b>enventualities </b>
</span> :

- ReverseMR : Reverse causality investigation
- PlotsMR : Graphical illustration of the results of an MR study
- add option to remove missense exonic variants (prots with aptamers)
