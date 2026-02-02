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
#> Package "ChIPseeker" needed for this function to work. Please install it.
#> Error in abs(cp.filtered.annotated$distanceToTSS): non-numeric argument to mathematical function

heatmap_differential(results)
#> Error: object 'results' not found
# }
```
