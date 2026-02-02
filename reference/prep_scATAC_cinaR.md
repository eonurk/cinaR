# prep_scATAC_cinaR

Prepare 10x scATAC peak-by-cell matrices for cinaR by pseudobulking per
sample.

## Usage

``` r
prep_scATAC_cinaR(
  counts,
  cell.meta,
  sample.col,
  group.col,
  cluster.col = NULL,
  peak.bed = NULL,
  min.cells = 20,
  verbose = TRUE
)
```

## Arguments

- counts:

  peak-by-cell count matrix (dense matrix or sparse dgCMatrix).

- cell.meta:

  data.frame with rownames as cell barcodes.

- sample.col:

  column name in cell.meta indicating biological replicate.

- group.col:

  column name in cell.meta indicating condition/group.

- cluster.col:

  optional column name for cell type/cluster. If provided, output is a
  named list per cluster.

- peak.bed:

  optional data.frame with CHR/START/STOP columns for peaks. If not
  provided, rownames(counts) are parsed as "chr:start-end" or
  "chr_start_end".

- min.cells:

  minimum number of cells required per sample (and per cluster if used).

- verbose:

  logical, prints informative messages.

## Value

list with elements \`bed\`, \`contrasts\`, and \`group.info\`, or a
named list of such lists when cluster.col is provided.

## Examples

``` r
# \donttest{
counts <- matrix(c(1, 0, 2, 1, 0, 1, 3, 0, 0, 2, 1, 0),
                 nrow = 2, byrow = TRUE)
rownames(counts) <- c("chr1:1-100", "chr1:101-200")
colnames(counts) <- paste0("cell", 1:6)
meta <- data.frame(sample = c("S1", "S1", "S2", "S2", "S3", "S3"),
                   group = c("A", "A", "B", "B", "B", "B"),
                   row.names = colnames(counts))
prep <- prep_scATAC_cinaR(counts, meta, sample.col = "sample", group.col = "group",
                          min.cells = 2)
# }
```
