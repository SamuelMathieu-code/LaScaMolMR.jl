<div align="center">

# LaScaMol MR

[![Build Status](https://github.com/SamuelMathieu-code/MrPainter.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/SamuelMathieu-code/MrPainter.jl/actions/workflows/CI.yml?query=branch%3Amain)

</div>

## Overview

LaScaMol.jl is a distrbuted Mendelian Randomization (MR) package that is focused on the generation of transcriptome wide / molecular MR analysies.

<span style="color:green">
<b>coded </b>
</span> :

- inputs : modular inputs of GWAS and QTL studies.
- ld : computation of the ld composite
- MrPref : implementation of heterogeneity tests and iv-regression methods.

<span style="color:yellow">
<b>ungoing </b>
</span> :

- MrStudyCis : Input parallel reading and filtering
- NaiveCis : MR study : clumping + MR

<span style="color:red">
<b>planned </b>
</span> :

- StrictCis : Strict approach to MR Study to handle local internal pleiotropy
- SecondChanceCis : alternative aproach to handle local internal pleiotropy

<span style="color:purple">
<b>enventualities </b>
</span> :

- TransMR : trans selection of ivs
- ReverseMR : Reverse causality investigation
- PlotsMR : Graphical illustration of the results of an MR study
- add option to remove missense exonic variants (prots with aptamers)

## TODO

- [ ] Implement & test NaiveCis.jl
    - [ ] Write unit test NaiveCis.jl
- [ ] Implement & test ld.jl
    - [x] test with real data in notebook and compare results with plink clumping
    - [ ] Write unit ld.jl
- [ ] Implement & test mrStudyCis.jl
    - [x] Implement harmonisation & indels removal
    - [ ] test with real data in notebook
    - [ ] Write unit test mrStudyCis.jl
- [ ] Benchmarks
    - [ ] Benchmark between @threads with collect(enumerate) and ThreadX.foreach (NaiveCis.jl)
