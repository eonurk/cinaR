#' prep_scATAC_cinaR
#'
#' Prepare 10x scATAC peak-by-cell matrices for cinaR by pseudobulking per sample.
#'
#' @param counts peak-by-cell count matrix (dense matrix or sparse dgCMatrix).
#' @param cell.meta data.frame with rownames as cell barcodes.
#' @param sample.col column name in cell.meta indicating biological replicate.
#' @param group.col column name in cell.meta indicating condition/group.
#' @param cluster.col optional column name for cell type/cluster. If provided,
#' output is a named list per cluster.
#' @param peak.bed optional data.frame with CHR/START/STOP columns for peaks.
#' If not provided, rownames(counts) are parsed as "chr:start-end" or "chr_start_end".
#' @param min.cells minimum number of cells required per sample (and per cluster if used).
#' @param verbose logical, prints informative messages.
#'
#' @return list with elements `bed`, `contrasts`, and `group.info`, or a named list
#' of such lists when cluster.col is provided.
#' @examples
#' \donttest{
#' counts <- matrix(c(1, 0, 2, 1, 0, 1, 3, 0, 0, 2, 1, 0),
#'                  nrow = 2, byrow = TRUE)
#' rownames(counts) <- c("chr1:1-100", "chr1:101-200")
#' colnames(counts) <- paste0("cell", 1:6)
#' meta <- data.frame(sample = c("S1", "S1", "S2", "S2", "S3", "S3"),
#'                    group = c("A", "A", "B", "B", "B", "B"),
#'                    row.names = colnames(counts))
#' prep <- prep_scATAC_cinaR(counts, meta, sample.col = "sample", group.col = "group",
#'                           min.cells = 2)
#' }
#' @export
prep_scATAC_cinaR <- function(counts,
                              cell.meta,
                              sample.col,
                              group.col,
                              cluster.col = NULL,
                              peak.bed = NULL,
                              min.cells = 20,
                              verbose = TRUE) {
  if (is.data.frame(counts)) {
    counts <- as.matrix(counts)
  }
  if (is.null(colnames(counts))) {
    stop("`counts` must have column names (cell barcodes).")
  }
  if (!is.data.frame(cell.meta)) {
    stop("`cell.meta` must be a data.frame with rownames as cell barcodes.")
  }
  if (is.null(rownames(cell.meta))) {
    stop("`cell.meta` must have rownames set to cell barcodes.")
  }
  if (!sample.col %in% colnames(cell.meta)) {
    stop("`sample.col` is not a column in `cell.meta`.")
  }
  if (!group.col %in% colnames(cell.meta)) {
    stop("`group.col` is not a column in `cell.meta`.")
  }

  common <- intersect(colnames(counts), rownames(cell.meta))
  if (length(common) == 0) {
    stop("No overlapping barcodes between `counts` and `cell.meta` rownames.")
  }
  if (length(common) < ncol(counts) && verbose) {
    message("Dropping ", ncol(counts) - length(common), " cells not found in `cell.meta`.")
  }
  counts <- counts[, common, drop = FALSE]
  cell.meta <- cell.meta[common, , drop = FALSE]

  if (!is.null(cluster.col)) {
    if (!cluster.col %in% colnames(cell.meta)) {
      stop("`cluster.col` is not a column in `cell.meta`.")
    }
    clusters <- unique(cell.meta[[cluster.col]])
    clusters <- clusters[!is.na(clusters)]
    res <- lapply(clusters, function(cl) {
      idx <- cell.meta[[cluster.col]] == cl
      .prep_sc_one(counts[, idx, drop = FALSE],
                   cell.meta[idx, , drop = FALSE],
                   sample.col,
                   group.col,
                   peak.bed,
                   min.cells,
                   verbose)
    })
    names(res) <- as.character(clusters)
    return(res)
  }

  .prep_sc_one(counts, cell.meta, sample.col, group.col, peak.bed, min.cells, verbose)
}

.prep_sc_one <- function(counts,
                         cell.meta,
                         sample.col,
                         group.col,
                         peak.bed,
                         min.cells,
                         verbose) {
  sample_ids <- as.character(cell.meta[[sample.col]])
  group_ids <- as.character(cell.meta[[group.col]])

  if (any(is.na(sample_ids)) || any(is.na(group_ids))) {
    keep <- !is.na(sample_ids) & !is.na(group_ids)
    if (sum(keep) == 0) {
      stop("All cells have NA in `sample.col` or `group.col`.")
    }
    if (verbose) {
      message("Dropping ", sum(!keep), " cells with NA sample/group.")
    }
    counts <- counts[, keep, drop = FALSE]
    sample_ids <- sample_ids[keep]
    group_ids <- group_ids[keep]
  }

  sample_counts <- table(sample_ids)
  keep_samples <- names(sample_counts)[sample_counts >= min.cells]
  if (length(keep_samples) == 0) {
    stop("No samples have at least `min.cells` cells.")
  }
  keep_cells <- sample_ids %in% keep_samples
  if (sum(!keep_cells) > 0 && verbose) {
    message("Dropping ", sum(!keep_cells), " cells from samples with < ", min.cells, " cells.")
  }
  counts <- counts[, keep_cells, drop = FALSE]
  sample_ids <- sample_ids[keep_cells]
  group_ids <- group_ids[keep_cells]

  sample_levels <- unique(sample_ids)

  group_by_sample <- tapply(group_ids, sample_ids, function(x) unique(x))
  bad_samples <- names(Filter(function(x) length(x) != 1, group_by_sample))
  if (length(bad_samples) > 0) {
    stop("Samples with multiple groups detected: ", paste(bad_samples, collapse = ", "))
  }

  pb <- .pseudobulk_counts(counts, sample_ids, sample_levels)

  group_info <- data.frame(
    sample = sample_levels,
    group = vapply(group_by_sample[sample_levels], function(x) x[1], character(1)),
    n_cells = as.integer(sample_counts[sample_levels]),
    row.names = sample_levels,
    stringsAsFactors = FALSE
  )

  if (length(unique(group_info$group)) < 2 && verbose) {
    message("Only one group detected after filtering; differential analysis may not be meaningful.")
  }

  peak_df <- .parse_peak_bed(rownames(counts), peak.bed)

  pb_mat <- as.matrix(pb)
  bed <- cbind(peak_df, pb_mat)
  contrasts <- group_info$group
  names(contrasts) <- rownames(group_info)

  list(bed = bed, contrasts = contrasts, group.info = group_info)
}

.pseudobulk_counts <- function(counts, sample_ids, sample_levels) {
  if (inherits(counts, "dgCMatrix") || inherits(counts, "Matrix")) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Matrix package is required for sparse inputs.")
    }
    group_fac <- factor(sample_ids, levels = sample_levels)
    mm <- Matrix::sparse.model.matrix(~0 + group_fac)
    colnames(mm) <- sub("^group_fac", "", colnames(mm))
    pb <- counts %*% mm
    colnames(pb) <- colnames(mm)
    return(pb)
  }

  pb <- vapply(sample_levels, function(s) {
    rowSums(counts[, sample_ids == s, drop = FALSE])
  }, numeric(nrow(counts)))
  rownames(pb) <- rownames(counts)
  pb
}

.parse_peak_bed <- function(peak_ids, peak.bed) {
  if (!is.null(peak.bed)) {
    if (inherits(peak.bed, "GRanges")) {
      peak.bed <- as.data.frame(peak.bed)
    }
    if (!is.data.frame(peak.bed)) {
      stop("`peak.bed` must be a data.frame or GRanges.")
    }
    if (nrow(peak.bed) != length(peak_ids)) {
      stop("`peak.bed` rows must match number of peaks in `counts`.")
    }
    cn <- tolower(colnames(peak.bed))
    chr_idx <- match("chr", cn)
    start_idx <- match("start", cn)
    end_idx <- match("end", cn)
    if (any(is.na(c(chr_idx, start_idx, end_idx)))) {
      peak_df <- peak.bed[, 1:3, drop = FALSE]
      colnames(peak_df) <- c("CHR", "START", "STOP")
      return(peak_df)
    }
    peak_df <- peak.bed[, c(chr_idx, start_idx, end_idx), drop = FALSE]
    colnames(peak_df) <- c("CHR", "START", "STOP")
    return(peak_df)
  }

  if (is.null(peak_ids)) {
    stop("Rownames of `counts` are required to parse peaks.")
  }

  example <- peak_ids[1]
  if (grepl(":", example, fixed = TRUE) && grepl("-", example, fixed = TRUE)) {
    parts <- strsplit(peak_ids, ":", fixed = TRUE)
    chr <- vapply(parts, `[`, character(1), 1)
    rest <- vapply(parts, `[`, character(1), 2)
    pos <- strsplit(rest, "-", fixed = TRUE)
    start <- vapply(pos, `[`, character(1), 1)
    end <- vapply(pos, `[`, character(1), 2)
  } else if (grepl("_", example, fixed = TRUE)) {
    parts <- strsplit(peak_ids, "_", fixed = TRUE)
    chr <- vapply(parts, `[`, character(1), 1)
    start <- vapply(parts, `[`, character(1), 2)
    end <- vapply(parts, `[`, character(1), 3)
  } else {
    stop("Peak IDs must be in 'chr:start-end' or 'chr_start_end' format.")
  }

  peak_df <- data.frame(
    CHR = chr,
    START = as.integer(start),
    STOP = as.integer(end),
    stringsAsFactors = FALSE
  )

  if (any(is.na(peak_df$START)) || any(is.na(peak_df$STOP))) {
    stop("Failed to parse peak coordinates from rownames.")
  }
  peak_df
}
