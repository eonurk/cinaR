# prep_scATAC_seurat

Prepare 10x scATAC data from a Seurat/Signac object for cinaR.

## Usage

``` r
prep_scATAC_seurat(
  object,
  assay = NULL,
  slot = "counts",
  sample.col,
  group.col,
  cluster.col = NULL,
  peak.bed = NULL,
  min.cells = 20,
  verbose = TRUE
)
```

## Arguments

- object:

  Seurat object containing an ATAC assay (typically "peaks").

- assay:

  assay name to use; defaults to Seurat::DefaultAssay(object).

- slot:

  assay slot to pull counts from (default "counts").

- sample.col:

  column name in object@meta.data indicating biological replicate.

- group.col:

  column name in object@meta.data indicating condition/group.

- cluster.col:

  optional column name for cell type/cluster.

- peak.bed:

  optional data.frame with CHR/START/STOP columns for peaks.

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
if (requireNamespace("Seurat", quietly = TRUE)) {
  prep <- prep_scATAC_seurat(seurat_obj,
                             sample.col = "sample",
                             group.col = "group",
                             assay = "peaks")
}
#> Error: object 'seurat_obj' not found
# }
```
