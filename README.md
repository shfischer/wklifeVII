length-based empirical data-limited catch rule
================

# FLR data-limited MSE for ICES WKLIFE

## Introduction

This repository contains the Management Strategy Evaluation (MSE) for
the ICES data-limited catch rule as presented during ICES WKLIFE VII,
VIII and IX. The simulation is based on the Fisheries Library in R
([FLR](http://www.flr-project.org/)) and makes use of the Assessment for
All (a4a) standard MSE framework ([`FLR/mse`](github.com/FLR/mse))
developed during the Workshop on development of MSE algorithms with
R/FLR/a4a ([Jardim et
al., 2017](https://ec.europa.eu/jrc/en/publication/assessment-all-initiativea4a-workshop-development-mse-algorithms-rflra4a)).

The repository contains the source code for the work published in:

> Simon H. Fischer, José A. A. De Oliveira, Laurence T. Kell (2020).
> Linking the performance of a data-limited empirical catch rule to
> life-history traits, ICES Journal of Marine Science,
> <https://doi.org/10.1093/icesjms/fsaa054>.

The state of the code for the publication is stored in release v1.0
(<https://github.com/shfischer/wklifeVII/releases/tag/v1.0>).

## Repository structure

The repository contains the following R scripts in the `R/` directory:

  - `OM1.R` & `OM2.R`: Scripts for creating the operating models for 29
    data-limited fish stocks,
  - `MP.R`: script for running the MSE scenarios and is called from a
    job submission script,
  - `MP_stats.R`: script for post processing the results from `MP.R`,
  - `MP_analysis.R`: script for analysing the results,
  - `MP_plots.R`: script for creating plots,
  - `MP_functions.R`: script with additional functions, used for
    creating the operaing models, run the MSE and processing it
    afterwards,
  - `input/`: contains csv files with the life-history parameters used
    to create the operating models

## R, R packages and version info

The MSE simulation was run on a high performance computing cluster:

``` r
sessionInfo()
R version 3.5.1 (2018-07-02) -- "Feather Spray"
Copyright (C) 2018 The R Foundation for Statistical Computing
Platform: x86_64-conda_cos6-linux-gnu (64-bit)
(...)
other attached packages:
 [1] doRNG_1.7.1        rngtools_1.3.1     pkgmaker_0.27      registry_0.5
 [5] mseDL_0.9.9        foreach_1.4.7      data.table_1.12.2  FLBRP_2.5.3
 [9] ggplotFL_2.6.6     ggplot2_3.2.1      FLAssess_2.6.3     FLash_2.5.11
[13] FLCore_2.6.11.9001 iterators_1.0.12   lattice_0.20-38
```

The framework uses FLR and requires the following FLR packages:

  - `FLCore`
  - `FLash`
  - `FLBRP`
  - `ggplotFL`
  - `FLife`
  - `mseDL` (a fork of the FLR/mse package for data-limited MSE)

The specific FLR package versions as used for the simulation can be
installed with `devtools`:

``` r
devtools::install_github(repo = "flr/FLCore", ref = "d55bc6570c0134c6bea6c3fc44be20378691e042")
devtools::install_github(repo = "flr/FLash", ref = "7c47560cf57627068259404bb553f2b644682726")
devtools::install_github(repo = "flr/FLBRP", ref = "142d5e14137c5ceb4526afd6718c26269ad81e7c")
devtools::install_github(repo = "flr/ggplotFL", ref = "9b502a1aa01524637f4f269a3353a92c7d452db0")
devtools::install_github(repo = "flr/FLife", ref = "d0cca5e574a77fb52ec607a25c244969b9c8dd38")
devtools::install_github(repo = "shfischer/mse", ref = "80b5cf18dc9611f7307f599564ccdfbad433948d")
```

Furthermore, some more R packages available from CRAN are required:

``` r
install.packages(c("foreach", "doParallel", "dplyr", "tidyr", "data.table",
                   "doRNG", "dtwclust", "ggdendro", "glmnet")) 
```
