# run_enrichment

This function is run, if the enrichment pipeline wants to be called
afterwards. Setting reference genome to the same genome which cinaR was
run should be given to this function!

## Usage

``` r
run_enrichment(
  results,
  geneset = NULL,
  experiment.type = "ATAC-Seq",
  reference.genome = NULL,
  enrichment.method = NULL,
  enrichment.FDR.cutoff = 1,
  background.genes.size = 20000,
  verbose = TRUE
)
```

## Arguments

- results:

  list, DA peaks list for different contrasts

- geneset:

  Pathways to be used in enrichment analyses. If not set vp2008
  (Chaussabel, 2008) immune modules will be used. This can be set to any
  geneset using \`read.gmt\` function from \`qusage\` package. Different
  modules are available: https://www.gsea-msigdb.org/gsea/downloads.jsp.

- experiment.type:

  The type of experiment either set to "ATAC-Seq" or "RNA-Seq"

- reference.genome:

  genome of interested species. It should be 'hg38', 'hg19' or 'mm10'.

- enrichment.method:

  There are two methodologies for enrichment analyses, Hyper-geometric
  p-value (HPEA) or Geneset Enrichment Analyses (GSEA).

- enrichment.FDR.cutoff:

  FDR cut-off for enriched terms, p-values are corrected by
  Benjamini-Hochberg procedure

- background.genes.size:

  number of background genes for hyper-geometric p-value calculations.
  Default is 20,000.

- verbose:

  prints messages through running the pipeline

## Value

list, enrichment analyses results along with corresponding differential
analyses outcomes

## Examples

``` r
# \donttest{
library(cinaR)
data(atac_seq_consensus_bm) # calls 'bed'

# a vector for comparing the examples
contrasts <- sapply(strsplit(colnames(bed), split = "-", fixed = TRUE),
                    function(x){x[1]})[4:25]

results <- cinaR(bed, contrasts, reference.genome = "mm10", run.enrichment = FALSE)
#> >> Experiment type: ATAC-Seq
#> >> Matrix is filtered!
#> >> preparing features information...        2026-02-02 12:28:11 
#> >> Using Genome: mm10 ...
#> >> identifying nearest features...      2026-02-02 12:28:11 
#> >> calculating distance from peak to TSS...     2026-02-02 12:28:11 
#> >> assigning genomic annotation...      2026-02-02 12:28:11 
#> >> Using Genome: mm10 ...
#> >> Using Genome: mm10 ...
#> >> assigning chromosome lengths             2026-02-02 12:28:13 
#> >> done...                  2026-02-02 12:28:13 
#> >> Method: edgeR
#>  FDR:0.05& abs(logFC)<0
#> >> Estimating dispersion...
#> >> Fitting GLM...
#> >> DA peaks are found!

results_with_enrichment <- run_enrichment(results, reference.genome = "mm10")
#> >> No `geneset` is specified so immune modules (Chaussabel, 2008) will be used!
#> >> enrichment.method` is not selected. Hyper-geometric p-value (HPEA) will be used!
#> >> Mice gene symbols are converted to human symbols!
# }
```
