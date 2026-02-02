#' prep_scATAC_seurat
#'
#' Prepare 10x scATAC data from a Seurat/Signac object for cinaR.
#'
#' @param object Seurat object containing an ATAC assay (typically "peaks").
#' @param assay assay name to use; defaults to Seurat::DefaultAssay(object).
#' @param slot assay slot to pull counts from (default "counts").
#' @param sample.col column name in object@meta.data indicating biological replicate.
#' @param group.col column name in object@meta.data indicating condition/group.
#' @param cluster.col optional column name for cell type/cluster.
#' @param peak.bed optional data.frame with CHR/START/STOP columns for peaks.
#' @param min.cells minimum number of cells required per sample (and per cluster if used).
#' @param verbose logical, prints informative messages.
#'
#' @return list with elements `bed`, `contrasts`, and `group.info`, or a named list
#' of such lists when cluster.col is provided.
#' @examples
#' \donttest{
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   prep <- prep_scATAC_seurat(seurat_obj,
#'                              sample.col = "sample",
#'                              group.col = "group",
#'                              assay = "peaks")
#' }
#' }
#' @export
prep_scATAC_seurat <- function(object,
                               assay = NULL,
                               slot = "counts",
                               sample.col,
                               group.col,
                               cluster.col = NULL,
                               peak.bed = NULL,
                               min.cells = 20,
                               verbose = TRUE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package \"Seurat\" needed for this function to work. Please install it.")
  }
  if (is.null(object)) {
    stop("`object` must be a Seurat object.")
  }
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }

  counts <- Seurat::GetAssayData(object, assay = assay, slot = slot)
  meta <- object@meta.data
  if (is.null(meta)) {
    stop("`object` does not contain meta.data.")
  }

  if (is.null(peak.bed)) {
    assay_obj <- object[[assay]]
    if (!is.null(assay_obj) && "ranges" %in% methods::slotNames(assay_obj)) {
      gr <- methods::slot(assay_obj, "ranges")
      if (inherits(gr, "GRanges")) {
        peak.bed <- data.frame(
          CHR = as.character(GenomicRanges::seqnames(gr)),
          START = GenomicRanges::start(gr),
          STOP = GenomicRanges::end(gr),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  prep_scATAC_cinaR(
    counts = counts,
    cell.meta = meta,
    sample.col = sample.col,
    group.col = group.col,
    cluster.col = cluster.col,
    peak.bed = peak.bed,
    min.cells = min.cells,
    verbose = verbose
  )
}
