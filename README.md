
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cinaR <a href='https://eonurk.github.io/cinaR'><img src='man/figures/cinaR.png' align="right" height="139" /></a>

<!-- badges: start -->

[![Build
Status](https://travis-ci.com/eonurk/cinaR.svg?branch=master)](https://travis-ci.com/eonurk/cinaR)
[![R build
status](https://github.com/eonurk/cinaR/workflows/R-CMD-check/badge.svg)](https://github.com/eonurk/cinaR/actions)
<!-- badges: end -->

## Overview

`cinaR` is a single wrapper function for end-to-end computational
analyses of bulk ATAC-seq profiles. It starts from a consensus peak file
and outputs Differentially Accessible (DA) peaks and Enrichment Analyses
results.

## Installation

``` r
library(devtools)
install_github("eonurk/cinaR")
```

## Quick Start

``` r
library(cinaR)
data("atac_seq_consensus_bm")

contrasts<- c("B6", "B6", "B6", "B6", "B6", "NZO", "NZO", "NZO", "NZO", "NZO", "NZO", 
              "B6", "B6", "B6", "B6", "B6", "NZO", "NZO", "NZO", "NZO", "NZO", "NZO")

results <- cinaR(bed, contrasts, reference.genome = "mm10")
#> 
#> >> preparing features information...      2020-10-19 16:42:44 
#> >> identifying nearest features...        2020-10-19 16:42:46 
#> >> calculating distance from peak to TSS...   2020-10-19 16:42:47 
#> >> assigning genomic annotation...        2020-10-19 16:42:47 
#> >> assigning chromosome lengths           2020-10-19 16:43:19 
#> >> done...                    2020-10-19 16:43:19 
#> >> Method: edgeR
#>  FDR: 0.05 & abs(logFC)< 0 
#> >> Estimating dispersion...
#> >> Fitting GLM...
#> >> DA peaks are found!
#> >> No `geneset` is specified so immune modules (Chaussabel, 2008) will be used!
#> >> enrichment.method` is not selected. Hyper-geometric p-value (HPEA) will be used!
#> >> Converting mouse genes to human...
#> >> Mouse to human mapping is created for all genes!
#> >> Human gene symbols are converted to mice!
#> >> Enrichment results are ready...
#> >> Done!
```

``` r
pca_plot(results, contrasts, show.names = F)
```

<img src="man/figures/unnamed-chunk-4-1.png" width="100%" />

``` r
dot_plot(results)
#> Warning: Removed 54 rows containing missing values (geom_point).
```

<img src="man/figures/unnamed-chunk-5-1.png" width="100%" />
