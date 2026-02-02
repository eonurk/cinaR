# heatmap_differential

plot differentially accessible peaks for a given comparison

## Usage

``` r
heatmap_differential(results, comparison = NULL, ...)
```

## Arguments

- results:

  cinaR result object

- comparison:

  these are created by cinaR from \`contrasts\` user provided. If not
  selected the first comparison will be shown!

- ...:

  additional arguments for heatmap function, for more info \`?pheatmap\`

## Value

ggplot object

## Examples

``` r
# \donttest{
library(cinaR)
data(atac_seq_consensus_bm) # calls 'bed'

# a vector for comparing the examples
contrasts <- sapply(strsplit(colnames(bed), split = "-", fixed = TRUE),
                    function(x){x[1]})[4:25]

results <- cinaR(bed, contrasts, reference.genome = "mm10")
#> >> Experiment type: ATAC-Seq
#> >> Matrix is filtered!
#> >> preparing features information...        2026-02-02 12:19:43 
#> >> Using Genome: mm10 ...
#> >> identifying nearest features...      2026-02-02 12:19:43 
#> >> calculating distance from peak to TSS...     2026-02-02 12:19:43 
#> >> assigning genomic annotation...      2026-02-02 12:19:43 
#> >> Using Genome: mm10 ...
#> >> Using Genome: mm10 ...
#> >> assigning chromosome lengths             2026-02-02 12:19:46 
#> >> done...                  2026-02-02 12:19:46 
#> >> Method: edgeR
#>  FDR:0.05& abs(logFC)<0
#> >> Estimating dispersion...
#> >> Fitting GLM...
#> >> DA peaks are found!
#> >> No `geneset` is specified so immune modules (Chaussabel, 2008) will be used!
#> >> enrichment.method` is not selected. Hyper-geometric p-value (HPEA) will be used!
#> >> Mice gene symbols are converted to human symbols!
#> >> Enrichment results are ready...
#> >> Done!

heatmap_differential(results)
#> Warning: 'comparison' is not set so 'B6_NZO' will be used!

# }
```
