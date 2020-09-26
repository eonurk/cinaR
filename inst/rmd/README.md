
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cinaR

<img align="left" width="125" height="125" src="../logo/cinaR.png">

   
`cinaR` is a single wrapper function for end-to-end computational
analyses of bulk ATAC-seq profiles. It starts from a consensus peak file
and outputs Differentially Accessible (DA) peaks with corresponding
genenames, adjusted p-values, log-fold changes and so on.    

## Installation

``` r
library(devtools)
install_github("eonurk/cinaR")
```

## Quick Start

``` r
library(cinaR)
data("atac_seq_consensus_bm")
```

Bed formatted consensus matrix (chr, start, end and samples)

``` r
dim(bed)
## [1] 1000   25
```

``` r
# bed formatted file
head(bed[,1:4])
##         Chr     Start       End B6-18mo-M-BM-47-GT18-01783
## 52834  chr5  24841478  24845196                       1592
## 29780 chr17   8162955   8164380                        109
## 67290  chr8  40577584  40578029                         72
## 51295  chr4 145277698 145278483                        110
## 4267   chr1 180808752 180815472                       2452
## 45102  chr3  88732151  88732652                         49
```

Create the contrasts you want to compare, here we create contrasts for
22 mice samples from different strains.

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
## >> preparing features information...      2020-09-11 11:10:18 
## >> identifying nearest features...        2020-09-11 11:10:19 
## >> calculating distance from peak to TSS...   2020-09-11 11:10:21 
## >> assigning genomic annotation...        2020-09-11 11:10:21 
## >> assigning chromosome lengths           2020-09-11 11:10:45 
## >> done...                    2020-09-11 11:10:45 
## >> Estimating dispersion...
## >> Fitting GLM...
## >> Method: edgeR
##  FDR: 0.05 & abs(logFC)< 0 
## >> DA peaks are found!
```

There are many information `cinaR` provides such as adjusted p value,
log fold-changes, gene names etc.

``` r
# B6 vs NZO
colnames(results$B6_NZO)
##  [1] "Row.names"     "seqnames"      "start"         "end"          
##  [5] "width"         "strand"        "annotation"    "geneChr"      
##  [9] "geneStart"     "geneEnd"       "geneLength"    "geneStrand"   
## [13] "geneId"        "transcriptId"  "distanceToTSS" "gene_name"    
## [17] "logFC"         "logCPM"        "F"             "PValue"       
## [21] "FDR"
```

Here is an overview of the result:

``` r
head(results$B6_NZO[,1:5])
##                  Row.names seqnames     start       end width
## 1 chr1_104737782_104738294     chr1 104737782 104738294   513
## 2 chr1_104738589_104739011     chr1 104738589 104739011   423
## 3 chr1_105662857_105665310     chr1 105662857 105665310  2454
## 4 chr1_105989457_105992240     chr1 105989457 105992240  2784
## 5 chr1_106153994_106154441     chr1 106153994 106154441   448
## 6 chr1_106170912_106173892     chr1 106170912 106173892  2981
```

## Creating different contrasts

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

## Selecting different reference genomes

For now, `cinaR` supports 3 genomes,for human and mice models:

  - `hg39`
  - `hg38`
  - `mm10`

You can set your it using `reference.genome` argument.

## Batch Effect Correction

If you suspect your data have unknown batch effects, you can use:

``` r
cinaR(..., batch.correction = T)
```

This option will run [Surrogate Variable
Analysis](http://bioconductor.org/packages/release/bioc/html/sva.html)
(SVA) and try to adjust your data for unknown batch effects. If however,
you already know the batches of the samples, you can simply set the
`batch.information` argument as well:

``` r
# batch information should be number a vector where
# the length of it equals to the number of samples.
cinaR(..., batch.correction = T, batch.information = c(rep(0, 11), rep(1,11)))
```

> Reminder - In our example data we have 22 samples

## Saving DA peaks to excel

Setting `save.DA.peaks = TRUE` in `cinaR` function will create a
`DApeaks.xlsx` file in the current directory. This file includes all the
comparisons in different tabs. Additionally, you can set the path/name
of the file using `DA.peaks.path` argument after setting `save.DA.peaks
= TRUE`.

For instance,

``` r
results <- cinaR(bed, contrasts, reference.genome = "mm10", 
                 save.DA.peaks = T, DA.peaks.path = "./Peaks_mice.xlsx")
```

will create an excel file with name `Peaks_mice.xlsx` in the current
directory.

## Using different GLM algorithms

Currently, `cinaR` supports 4 different algorithms, namely;

  - edgeR
  - limma-voom
  - limma-trend
  - DESeq2

If not set, it uses `edgeR` for differential analyses. You can change
the used algorithm by simply setting `DA.choice` argument. For more
information, `?cinaR`

## Contribution

You can send pull requests to make your contributions.

I occasionally mess up, so all comments are appreciated\!

## Future work

  - Add enrichment pipeline
      - hyper-geometric p-value
      - geneset enrichment analyses
      - make it compitable with `.gmt` format

## Author

  - [E Onur Karakaslar](https://eonurk.github.io)

## License

  - GNU General Public License v3.0

## References

  - Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor
    package for differential expression analysis of digital gene
    expression data.” Bioinformatics, 26(1), 139-140. doi:
    10.1093/bioinformatics/btp616.

  - Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015).
    “limma powers differential expression analyses for RNA-sequencing
    and microarray studies.” Nucleic Acids Research, 43(7), e47.

  - Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of
    fold change and dispersion for RNA-seq data with DESeq2. Genome
    Biology, 15:550. 10.1186/s13059-014-0550-8
