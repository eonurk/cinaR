
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Quick Start

``` r
library(cinaR)
data("atac_seq_consensus_bm")
```

Bed formatted consensus matrix (chr, start, end and samples)

``` r
dim(bed)
## [1] 75281    25
```

``` r
# bed formatted file
head(bed[,1:5])
##    Chr   Start     End B6-18mo-M-BM-47-GT18-01783 B6-18mo-F-BM-42-GT18-01780
## 1 chr1 3819184 3819579                         15                         24
## 2 chr1 4671555 4672020                         56                         39
## 3 chr1 4770013 4770279                         31                         10
## 4 chr1 4774906 4775814                        133                        107
## 5 chr1 4780053 4780656                         66                         48
## 6 chr1 4785296 4786419                        494                        342
```

Create the contrasts you want to compare

``` r
# create contrast vector which will be compared.
contrasts <- sapply(strsplit(colnames(bed), split = "-", fixed = T), 
                    function(x){x[1]})[4:25]
contrasts
##  [1] "B6"  "B6"  "B6"  "B6"  "B6"  "NZO" "NZO" "NZO" "NZO" "NZO" "NZO" "B6" 
## [13] "B6"  "B6"  "B6"  "B6"  "NZO" "NZO" "NZO" "NZO" "NZO" "NZO"
```

`cinaR` function directly computes the differentially accessible peaks.

``` r
# If reference genome is not set hg38 will be used!
results <- cinaR(bed, contrasts, reference.genome = "mm10")
## 
## Registered S3 method overwritten by 'enrichplot':
##   method               from
##   fortify.enrichResult DOSE
## >> preparing features information...      2020-09-11 10:27:13 
## >> identifying nearest features...        2020-09-11 10:27:14 
## >> calculating distance from peak to TSS...   2020-09-11 10:27:17 
## >> assigning genomic annotation...        2020-09-11 10:27:17 
## >> assigning chromosome lengths           2020-09-11 10:27:40 
## >> done...                    2020-09-11 10:27:40 
## >> Estimating dispersion...
## >> Fitting GLM...
## >> Method: edgeR
##  FDR: 0.05 & abs(logFC)< 0 
## >> DA peaks are found!
```

``` r
colnames(results[[1]])
##  [1] "Row.names"     "seqnames"      "start"         "end"          
##  [5] "width"         "strand"        "annotation"    "geneChr"      
##  [9] "geneStart"     "geneEnd"       "geneLength"    "geneStrand"   
## [13] "geneId"        "transcriptId"  "distanceToTSS" "gene_name"    
## [17] "logFC"         "logCPM"        "F"             "PValue"       
## [21] "FDR"
```

``` r
head(results[[1]][,1:5])
##                  Row.names seqnames     start       end width
## 1 chr1_104737782_104738294     chr1 104737782 104738294   513
## 2 chr1_104738589_104739011     chr1 104738589 104739011   423
## 3 chr1_105662857_105665310     chr1 105662857 105665310  2454
## 4 chr1_105989457_105992240     chr1 105989457 105992240  2784
## 5 chr1_106153994_106154441     chr1 106153994 106154441   448
## 6 chr1_106170912_106173892     chr1 106170912 106173892  2981
```

# Creating contrasts

Note that you can further divide the resolution of contrasts, for
instance this is also a valid vector

``` r
contrasts <- sapply(strsplit(colnames(bed), split = "-", fixed = T), 
                    function(x){paste(x[1:4], collapse = "-")})[4:25]
unique(contrasts)
## [1] "B6-18mo-M-BM"  "B6-18mo-F-BM"  "NZO-18mo-F-BM" "NZO-18mo-M-BM"
## [5] "B6-3mo-F-BM"   "B6-3mo-M-BM"   "NZO-3mo-F-BM"  "NZO-3mo-M-BM"
```

in this case, each of them will be compared to each other which will
result in 28 different comparisons.
