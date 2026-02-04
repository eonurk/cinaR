# Introduction to cinaR

    ## 

## Quick Start

``` r
library(cinaR)
data("atac_seq_consensus_bm")
```

Bed formatted consensus matrix (chr, start, end and samples)

``` r
dim(bed)
```

    ## [1] 1000   25

``` r
# bed formatted file
head(bed[,1:4])
```

    ##         Chr     Start       End B6-18mo-M-BM-47-GT18-01783
    ## 52834  chr5  24841478  24845196                       1592
    ## 29780 chr17   8162955   8164380                        109
    ## 67290  chr8  40577584  40578029                         72
    ## 51295  chr4 145277698 145278483                        110
    ## 4267   chr1 180808752 180815472                       2452
    ## 45102  chr3  88732151  88732652                         49

Create the contrasts you want to compare, here we create contrasts for
22 mice samples from different strains.

``` r
# create contrast vector which will be compared.
contrasts<- c("B6", "B6", "B6", "B6", "B6", "NZO", "NZO", "NZO", "NZO", "NZO", "NZO", 
              "B6", "B6", "B6", "B6", "B6", "NZO", "NZO", "NZO", "NZO", "NZO", "NZO")
```

`cinaR` function directly computes the differentially accessible peaks.

``` r
# If reference genome is not set hg38 will be used!
results <- cinaR(bed, contrasts, reference.genome = "mm10")
```

Now, you can access differential accessibility (DA) and enrichment
results.

``` r
names(results)
```

Inside `DA.results`, you have the consensus peaks (cp) and
differentially accessible (DA) peaks. If batch correction was run, then
`cp` will be a batch-corrected consensus matrix, otherwise it is the
filtered and normalized version of original consensus peaks you
provided.

``` r
names(results$DA.results)
```

There are many information `cinaR` provides such as adjusted p value,
log fold-changes, gene names etc for each peak:

``` r
colnames(results$DA.results$DA.peaks$B6_NZO)
```

Here is an overview of those DA peaks:

``` r
head(results$DA.results$DA.peaks$B6_NZO[,1:5])
```

> Since the comparison is `B6_NZO`, if fold-changes are positive it
> means they are more accesible in B6 compared to NZO and vice versa for
> negative values!

and here is a little overview for enrichment analyses results:

``` r
head(results$Enrichment.Results$B6_NZO[,c("module.name","overlapping.genes", "adj.p")])
```

### PCA Plots

You can easily get the PCA plots of the samples:

``` r
pca_plot(results, contrasts, show.names = F)
```

You can overlay different information onto PCA plots as well!

``` r
# Overlaid information
overlaid.info <- c("B6-18mo", "B6-18mo", "B6-18mo", "B6-18mo", "B6-18mo", 
                   "NZO-18mo", "NZO-18mo", "NZO-18mo", "NZO-18mo", "NZO-18mo", "NZO-18mo", 
                   "B6-3mo", "B6-3mo", "B6-3mo", "B6-3mo", "B6-3mo", 
                   "NZO-3mo", "NZO-3mo", "NZO-3mo", "NZO-3mo", "NZO-3mo", "NZO-3mo")
# Sample IDs
sample.names <- c("S01783", "S01780", "S01781", "S01778", "S01779", 
"S03804", "S03805", "S03806", "S03807", "S03808", 
"S03809", "S04678", "S04679", "S04680", "S04681", 
"S04682", "S10918", "S10916", "S10919", "S10921", 
"S10917", "S10920")
```

``` r
pca_plot(results, overlaid.info, sample.names)
```

### Heatmaps

#### Differential peaks

You can see the available comparisons using:

``` r
show_comparisons(results)
```

Then, plot the differential peaks for a selected contrast using:

``` r
heatmap_differential(results, comparison = "B6_NZO")
```

Also, you can configure your heatmaps using the additional arguments of
`pheatmap` function. For more information check out `?pheatmap`.

``` r
heatmap_differential(results, comparison = "B6_NZO", show_colnames = FALSE)
```

#### Most variable peaks

You can also plot most variable 100 peaks for all samples:

``` r
heatmap_var_peaks(results)
```

Plus, you can set the number of peaks to be used in these plots, and
again you can change the additional arguments of `pheatmap` function.
For more information check out `?pheatmap`.

``` r
heatmap_var_peaks(results, heatmap.peak.count = 200, cluster_cols = F)
```

### Enrichment Plots

You can plot your enrichment results using:

``` r
dot_plot(results)
```

if it gets too crowded, you can get rid of the irrelevant pathways as
well:

``` r
dot_plot(results, filter.pathways = T)
```

## Creating different contrasts

Note that you can further divide the resolution of contrasts, for
instance this is also a valid vector

``` r
contrasts <- sapply(strsplit(colnames(bed), split = "-", fixed = T), 
                    function(x){paste(x[1:4], collapse = ".")})[4:25]
unique(contrasts)
```

    ## [1] "B6.18mo.M.BM"  "B6.18mo.F.BM"  "NZO.18mo.F.BM" "NZO.18mo.M.BM"
    ## [5] "B6.3mo.F.BM"   "B6.3mo.M.BM"   "NZO.3mo.F.BM"  "NZO.3mo.M.BM"

in this case, each of them will be compared to each other which will
result in 28 different comparisons.

## Comparison scheme

As default, `cinaR` will use one vs one (OVA) scheme, comparing each
contrast to others one by one. However, this can be changed to one vs
all (OVA) which will compare each contrast to rest:

```
# one-vs-one (results in n choose k comparisons, default)
cinaR(..., comparison.scheme = "OVO")

# one-vs-all (results in n comparisons)
cinaR(..., comparison.scheme = "OVA")
```

## Running for bulk RNA-seq data

To run `cinaR` with RNA-seq experiments, just provide the count matrix,
and specify the experiment type:

``` r
cinaR(matrix = count.matrix, ..., experiment.type = "RNA-Seq")
```

Note that, `count.matrix` should be in the form of $g \times (1 + n)$
where $g$ is the number of genes and $n$ is the number of samples, and
plus one is for gene names.

> Note that currently `cinaR` can only handle gene symbols (e.g. FOSL2,
> FOXA) and ensembl stable IDs (e.g. ENSG00000010404) for both human and
> mice!

## Single-cell ATAC-seq (10x scATAC) preprocessing

`cinaR` is designed for bulk ATAC-seq. For 10x scATAC, first pseudobulk
your peak-by-cell matrix to preserve biological replicates, then run
`cinaR` on the resulting consensus matrix.

``` r
set.seed(1)
counts <- matrix(rpois(4 * 8, lambda = 2), nrow = 4)
rownames(counts) <- c("chr1:1-100", "chr1:101-200", "chr2:1-150", "chr2:151-300")
colnames(counts) <- paste0("cell", 1:8)

meta <- data.frame(
  sample = c("S1", "S1", "S2", "S2", "S3", "S3", "S4", "S4"),
  group = c("A", "A", "A", "A", "B", "B", "B", "B"),
  row.names = colnames(counts),
  stringsAsFactors = FALSE
)

prep <- prep_scATAC_cinaR(
  counts,
  meta,
  sample.col = "sample",
  group.col = "group",
  min.cells = 2
)

dim(prep$bed)
```

    ## [1] 4 7

``` r
prep$contrasts
```

    ##  S1  S2  S3  S4 
    ## "A" "A" "B" "B"

You can then pass `prep$bed` and `prep$contrasts` to
[`cinaR()`](https://eonurk.github.io/cinaR/reference/cinaR.md) as usual.

If you use a Seurat/Signac object, you can prepare the data directly:

``` r
prep <- prep_scATAC_seurat(seurat_obj,
                           sample.col = "sample",
                           group.col = "group",
                           assay = "peaks")
```

## Running enrichment with different dataset

You can run the enrichment analyses with a custom gene set:

``` r
cinaR(..., geneset = new_geneset)
```

#### `cinaRgenesets`

Easiest way to do this is to use
[cinaRgenesets](https://github.com/eonurk/cinaR-genesets) package. You
can select your gene set of interest and just plug it into your
pipeline.

#### MSigDB

You can also download different gene sets from [MSigDB
site](https://www.gsea-msigdb.org/gsea/downloads.jsp). Note that you
should use the human gene symbol versions.

> You can use `read.gmt` function from `qusage` package to read these
> genesets into your current environment.

#### Custom gene sets

A `geneset` must be a `.gmt` formatted symbol file.  
You can familiarize yourself with the format by checking out :

``` r
# default geneset to be used
data("VP2008")
```

> If you have gene and pathway names in a `data.frame`, you can use
> `split` function to create your own `.gmt` formatted gene sets
> e.g. `split(df$genes, df$pathways)`.

## Selecting different reference genomes

For now, `cinaR` supports 3 genomes for human and mice models:

- `hg38`
- `hg19`
- `mm10`

You can set your it using `reference.genome` argument.

## Batch Effect Correction

If you suspect your data have unknown batch effects, you can use:

``` r
cinaR(..., batch.correction = T)
```

This option will run [Surrogate Variable
Analysis](https://bioconductor.org/packages/release/bioc/html/sva.html)
(SVA) and try to adjust your data for unknown batch effects. If however,
you already know the batches of the samples, you can simply set the
`batch.information` argument as well. This will not run the SVA but just
add the batches to design matrix.

``` r
# runs SVA
cinaR(..., batch.correction = T)

# runs SVA with 2 surrogate variables
cinaR(..., batch.correction = T, sv.number = 2)

# adds only batch information to the design matrix! (does not run SVA)
# batch.information should be number a vector where
# the length of it equals to the number of samples.
cinaR(..., batch.correction = T, batch.information = c(rep(0, 11), rep(1,11)))
```

> Reminder - In our example data we have 22 samples

## Adding extra covariates

Sometimes, one might want to add additional covariates to adjust the
design matrix further, which affects the down-stream analyses. For
instance, ages or sexes of the samples could be additional covariates.
To account for those:

``` r
# Ages of the samples could be not in biological interests but should be accounted for!
cinaR(..., additional.covariates = c((18, 11), (3, 11)))

# More than one covariate for instance, sex and age
sex.info <- c("M", "F", "M", "F", "F", "F", "F", "F", "M", "M", "M", 
              "F", "F", "M", "M", "M", "F", "F", "M", "M", "F", "M")

age.info <- c((18, 11), (3, 11)
covs <- data.frame(Sex = sex.info, Age = age.info)

cinaR(..., additional.covariates = covs)
```

## Saving DA peaks to excel

Setting `save.DA.peaks = TRUE` in `cinaR` function will create a
`DApeaks.xlsx` file in the current directory. This file includes all the
comparisons in different tabs. Additionally, you can set the path/name
of the file using `DA.peaks.path` argument after setting
`save.DA.peaks = TRUE`.

For instance,

``` r
results <- cinaR(..., save.DA.peaks = T, DA.peaks.path = "./Peaks_mice.xlsx")
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
information,
[`?cinaR`](https://eonurk.github.io/cinaR/reference/cinaR.md)

## Some useful arguments

``` r
# new FDR threshold for DA peaks
results <- cinaR(..., DA.fdr.threshold = 0.1)

# filters out pathways
results <- cinaR(..., enrichment.FDR.cutoff = 0.1)

# does not run enrichment pipeline
results <- cinaR(..., run.enrichment = FALSE)

# creates the piechart from chIpSeeker package
results <- cinaR(..., show.annotation.pie = TRUE)

# change cut-off value for dot plots
dot_plot(..., fdr.cutoff = 0.05)
```

## References

- Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor
  package for differential expression analysis of digital gene
  expression data.” Bioinformatics, 26(1), 139-140. doi:
  10.1093/bioinformatics/btp616.

- Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015).
  “limma powers differential expression analyses for RNA-sequencing and
  microarray studies.” Nucleic Acids Research, 43(7), e47.

- Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold
  change and dispersion for RNA-seq data with DESeq2. Genome Biology,
  15:550. 10.1186/s13059-014-0550-8

## Session Info

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] cinaR_0.2.6
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] DBI_1.2.3                               
    ##   [2] bitops_1.0-9                            
    ##   [3] rlang_1.1.7                             
    ##   [4] magrittr_2.0.4                          
    ##   [5] DOSE_4.4.0                              
    ##   [6] otel_0.2.0                              
    ##   [7] matrixStats_1.5.0                       
    ##   [8] compiler_4.5.2                          
    ##   [9] RSQLite_2.4.5                           
    ##  [10] GenomicFeatures_1.62.0                  
    ##  [11] png_0.1-8                               
    ##  [12] systemfonts_1.3.1                       
    ##  [13] vctrs_0.7.1                             
    ##  [14] reshape2_1.4.5                          
    ##  [15] stringr_1.6.0                           
    ##  [16] pkgconfig_2.0.3                         
    ##  [17] crayon_1.5.3                            
    ##  [18] fastmap_1.2.0                           
    ##  [19] XVector_0.50.0                          
    ##  [20] Rsamtools_2.26.0                        
    ##  [21] rmarkdown_2.30                          
    ##  [22] UCSC.utils_1.6.1                        
    ##  [23] ragg_1.5.0                              
    ##  [24] bit_4.6.0                               
    ##  [25] xfun_0.56                               
    ##  [26] cachem_1.1.0                            
    ##  [27] cigarillo_1.0.0                         
    ##  [28] aplot_0.2.9                             
    ##  [29] GenomeInfoDb_1.46.2                     
    ##  [30] jsonlite_2.0.0                          
    ##  [31] blob_1.3.0                              
    ##  [32] DelayedArray_0.36.0                     
    ##  [33] BiocParallel_1.44.0                     
    ##  [34] parallel_4.5.2                          
    ##  [35] R6_2.6.1                                
    ##  [36] stringi_1.8.7                           
    ##  [37] bslib_0.10.0                            
    ##  [38] RColorBrewer_1.1-3                      
    ##  [39] limma_3.66.0                            
    ##  [40] rtracklayer_1.70.1                      
    ##  [41] boot_1.3-32                             
    ##  [42] GenomicRanges_1.62.1                    
    ##  [43] jquerylib_0.1.4                         
    ##  [44] GOSemSim_2.36.0                         
    ##  [45] Rcpp_1.1.1                              
    ##  [46] Seqinfo_1.0.0                           
    ##  [47] SummarizedExperiment_1.40.0             
    ##  [48] knitr_1.51                              
    ##  [49] R.utils_2.13.0                          
    ##  [50] IRanges_2.44.0                          
    ##  [51] splines_4.5.2                           
    ##  [52] Matrix_1.7-4                            
    ##  [53] tidyselect_1.2.1                        
    ##  [54] qvalue_2.42.0                           
    ##  [55] abind_1.4-8                             
    ##  [56] yaml_2.3.12                             
    ##  [57] codetools_0.2-20                        
    ##  [58] curl_7.0.0                              
    ##  [59] plyr_1.8.9                              
    ##  [60] lattice_0.22-7                          
    ##  [61] tibble_3.3.1                            
    ##  [62] Biobase_2.70.0                          
    ##  [63] KEGGREST_1.50.0                         
    ##  [64] S7_0.2.1                                
    ##  [65] evaluate_1.0.5                          
    ##  [66] gridGraphics_0.5-1                      
    ##  [67] desc_1.4.3                              
    ##  [68] Biostrings_2.78.0                       
    ##  [69] BiocManager_1.30.27                     
    ##  [70] pillar_1.11.1                           
    ##  [71] MatrixGenerics_1.22.0                   
    ##  [72] TxDb.Hsapiens.UCSC.hg19.knownGene_3.22.1
    ##  [73] stats4_4.5.2                            
    ##  [74] ggfun_0.2.0                             
    ##  [75] generics_0.1.4                          
    ##  [76] RCurl_1.98-1.17                         
    ##  [77] S4Vectors_0.48.0                        
    ##  [78] ggplot2_4.0.2                           
    ##  [79] scales_1.4.0                            
    ##  [80] glue_1.8.0                              
    ##  [81] tools_4.5.2                             
    ##  [82] BiocIO_1.20.0                           
    ##  [83] ggnewscale_0.5.2                        
    ##  [84] data.table_1.18.2.1                     
    ##  [85] locfit_1.5-9.12                         
    ##  [86] fgsea_1.36.2                            
    ##  [87] GenomicAlignments_1.46.0                
    ##  [88] fs_1.6.6                                
    ##  [89] XML_3.99-0.20                           
    ##  [90] fastmatch_1.1-8                         
    ##  [91] cowplot_1.2.0                           
    ##  [92] grid_4.5.2                              
    ##  [93] edgeR_4.8.2                             
    ##  [94] AnnotationDbi_1.72.0                    
    ##  [95] patchwork_1.3.2                         
    ##  [96] restfulr_0.0.16                         
    ##  [97] cli_3.6.5                               
    ##  [98] rappdirs_0.3.4                          
    ##  [99] textshaping_1.0.4                       
    ## [100] S4Arrays_1.10.1                         
    ## [101] dplyr_1.2.0                             
    ## [102] gtable_0.3.6                            
    ## [103] R.methodsS3_1.8.2                       
    ## [104] yulab.utils_0.2.3                       
    ## [105] sass_0.4.10                             
    ## [106] digest_0.6.39                           
    ## [107] BiocGenerics_0.56.0                     
    ## [108] ggrepel_0.9.6                           
    ## [109] SparseArray_1.10.8                      
    ## [110] ggplotify_0.1.3                         
    ## [111] rjson_0.2.23                            
    ## [112] htmlwidgets_1.6.4                       
    ## [113] farver_2.1.2                            
    ## [114] memoise_2.0.1                           
    ## [115] htmltools_0.5.9                         
    ## [116] pkgdown_2.2.0                           
    ## [117] R.oo_1.27.1                             
    ## [118] lifecycle_1.0.5                         
    ## [119] httr_1.4.7                              
    ## [120] statmod_1.5.1                           
    ## [121] GO.db_3.22.0                            
    ## [122] bit64_4.6.0-1
