# dot_plot

Given the results from \`cinaR\` it produces dot plots for enrichment
analyses.

## Usage

``` r
dot_plot(results, fdr.cutoff = 0.1, filter.pathways = FALSE)
```

## Arguments

- results:

  cinaR result object

- fdr.cutoff:

  Pathways with smaller fdr values than the cut-off will be shown as
  dots.

- filter.pathways:

  logical, it will filter the pathways from dot plot with fdr values
  less than \`fdr.cutoff\`.

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
#> Package "ChIPseeker" needed for this function to work. Please install it.
#> Error in abs(cp.filtered.annotated$distanceToTSS): non-numeric argument to mathematical function

dot_plot(results)
#> Error: object 'results' not found
# }
```
