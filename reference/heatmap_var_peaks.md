# heatmap_var_peaks

plot most variable k peaks (default k = 100) among all samples

## Usage

``` r
heatmap_var_peaks(results, heatmap.peak.count = 100, ...)
```

## Arguments

- results:

  cinaR result object

- heatmap.peak.count:

  number of peaks to be plotted. If number of peaks are less than k then
  all peaks will be used.

- ...:

  additional arguments for heatmap function, for more info \`?pheatmap\`

## Value

ggplot object

## Examples

``` r
library(cinaR)
data(atac_seq_consensus_bm) # calls 'bed'

# creating dummy results
results <- NULL
results[["cp"]] <- bed[,c(4:25)]

heatmap_var_peaks(results)

```
