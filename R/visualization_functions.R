#' dot_plot
#'
#' Given the results from `cinaR` it produces dot plots for enrichment analyses.
#'
#' @param results cinaR result object
#' @param fdr.cutoff Pathways with smaller fdr values than the cut-off
#' will be shown as dots.
#' @param filter.pathways logical, it will filter the pathways from dot plot
#' with fdr values less than `fdr.cutoff`.
#' @return ggplot object
#'
#' @examples
#' \donttest{
#' library(cinaR)
#' data(atac_seq_consensus_bm) # calls 'bed'
#'
#' # a vector for comparing the examples
#' contrasts <- sapply(strsplit(colnames(bed), split = "-", fixed = TRUE),
#'                     function(x){x[1]})[4:25]
#'
#' results <- cinaR(bed, contrasts, reference.genome = "mm10")
#'
#' dot_plot(results)
#' }
#' @export
dot_plot <- function(results, fdr.cutoff = 0.1, filter.pathways = FALSE){

  if(is.null (results[["Enrichment.Results"]])){
    stop("Did you run the enrichment pipeline in cinaR? For more info ?cinaR")
  }

  #silence build NOTES
  adj.p <- column_label <-
    log2FoldChange <- module.name <- padj <-
    prcomp <- status <- NULL

  # Set ATAC-Seq and check if it is correct
  exp.atac <- TRUE

  atac.check <- unlist(lapply(results[["DA.results"]][["DA.peaks"]],
                       function(x){ncol(x) == 3}), use.names = FALSE)

  if(any(atac.check)){
    message(">> cinaR was run for RNA-Seq experiments!")
    exp.atac <- FALSE
  }

  results <- results[["Enrichment.Results"]]

  # add list names as a column
  df.plot <- dplyr::bind_rows(results, .id = "column_label")

  # GSEA was run on the data
  if("padj" %in% colnames(df.plot)){

    # create a status column for plot
    df.plot$status <- ifelse(df.plot$NES > 0, "Opening", "Closing")

    # change the required field names
    colnames(df.plot)[c(2,4)] <- c("module.name", "adj.p")
  }

  if (!exp.atac){
    df.plot[,"status"][ df.plot[,"status"] == "Opening"] <- "Up"
    df.plot[,"status"][ df.plot[,"status"] == "Closing"] <- "Down"
  }

  if(filter.pathways){
    if (sum(df.plot$adj.p < fdr.cutoff) == 0){
      stop("You can't filter because there are no pathways to be displayed!")
    }
    df.plot <- subset(df.plot, adj.p < fdr.cutoff)
  }

  # create ggplot
  plot.dot <- ggplot2::ggplot(df.plot,
                              ggplot2::aes(x = column_label,
                                           y = module.name,
                                           size = ifelse(adj.p < fdr.cutoff, -log(adj.p), NA),
                                           color = status))

  plot.dot <- plot.dot + ggplot2::geom_point()

  plot.dot <- plot.dot + ggplot2::labs(x = "Contrast",
                                       y = "Pathways",
                                       color = "Sign",
                                       size = "-log10(adj.p)",
                                       caption = paste0("FDR < ", fdr.cutoff))

  plot.dot <- plot.dot + ggplot2::scale_color_manual(values = color_values)

  plot.dot <- plot.dot + ggplot2::theme_minimal()

  return(plot.dot)
}


#' pca_plot
#'
#' @param results cinaR result object
#' @param overlaid.info overlaid information onto the samples
#' @param sample.names names of the samples shown on pca plot
#' @param show.names logical, if set FALSE sample names will be hidden
#' @return ggplot object
#'
#' @examples
#'
#' #' library(cinaR)
#' data(atac_seq_consensus_bm) # calls 'bed'
#'
#' # creating dummy results
#' results <- NULL
#' results[["cp"]] <- bed[,c(4:25)]
#'
#' # a vector for comparing the examples
#' contrasts <- sapply(strsplit(colnames(bed), split = "-", fixed = TRUE),
#'                     function(x){x[1]})[4:25]
#'
#' ## overlays the contrasts info onto PCA plots
#' pca_plot(results, contrasts)
#'
#' ## you can overlay other information as well,
#' ## as long as it is the same length with the
#' ## number of samples.
#'
#' sample.info <- c(rep("Group A", 11), rep("Group B", 11))
#' pca_plot(results, sample.info, show.names = FALSE)
#'
#' @export
pca_plot <- function(results, overlaid.info, sample.names = NULL, show.names = TRUE){

  #silence CRAN build NOTES
  PC1 <- PC2 <- NULL

  # if enrichment not run!
  if(!is.null(results[["cp"]])){
    cp <- results[["cp"]]
  } else { # if run
    cp <- results[["DA.results"]][["cp"]]
  }

  if(is.null(sample.names)){
    sample.names <- colnames(cp)
  } else{
    if(length(sample.names) != ncol(cp)){
      stop("The length of `sample.names` should be equal to number of samples.")
    }
  }

  # eliminate NaN values before-hand if there is any.
  pca <- stats::prcomp(t(stats::na.omit(cp)), center = TRUE)

  d  <- round(pca$sdev^2/sum(pca$sdev^2)*100, digits=1)
  xl <- sprintf("PC 1: %.1f %%", d[1])
  yl <- sprintf("PC 2: %.1f %%", d[2])


  plot.df <- data.frame(PC1 = as.numeric(pca$x[,1]),
                   PC2 = as.numeric(pca$x[,2]),
                   overlaid.info = overlaid.info,
                   names = sample.names
                   )

  plot.pca <- ggplot2::ggplot(plot.df, ggplot2::aes(PC1, PC2, color = overlaid.info)) +
    ggplot2::geom_point(size = 4) +
    ggplot2::labs(x=xl,y=yl) +
    ggplot2::theme_minimal() +
    ggplot2::labs(color = "Status") +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::theme_light()

  if (typeof(overlaid.info) %in% c("character", "factor")){
    plot.pca <- plot.pca +
      ggplot2::scale_color_manual(values = RColorBrewer::brewer.pal(n = 9, name = "Set1"))
  }

  if(show.names){
    plot.pca <- plot.pca + ggrepel::geom_text_repel(ggplot2::aes(label = names))
  }

  return(plot.pca)
}

#' heatmap_var_peaks
#'
#' plot most variable k peaks (default k = 100) among all samples
#'
#' @param results cinaR result object
#' @param heatmap.peak.count number of peaks to be plotted.
#' If number of peaks are less than k then all peaks will be used.
#' @param ... additional arguments for heatmap function, for more info `?pheatmap`
#' @return ggplot object
#' @examples
#' library(cinaR)
#' data(atac_seq_consensus_bm) # calls 'bed'
#'
#' # creating dummy results
#' results <- NULL
#' results[["cp"]] <- bed[,c(4:25)]
#'
#' heatmap_var_peaks(results)
#'
#' @export
heatmap_var_peaks <- function(results, heatmap.peak.count = 100, ...){

  # if enrichment not run!
  if(!is.null(results[["cp"]])){
    cp <- results[["cp"]]
  } else { # if run
    cp <- results[["DA.results"]][["cp"]]
  }

  cp <- stats::na.omit(cp)

  # Remove possible na's from data and
  # reorder peak according to their standard deviation in decreasing order
  cp <- cp [rev(order(apply(cp, 1, stats::sd))),]

  # normalize peak distributions
  # utils function
  mat.heatmap <- scale_rows(cp)

  mat.heatmap <- mat.heatmap[1:min(nrow(cp),heatmap.peak.count),]

  breaksList = seq(min(mat.heatmap), max(mat.heatmap), by = .01)
  plot.pheatmap <- pheatmap::pheatmap(mat.heatmap, scale = "none",
                     color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n=7,"RdBu")))(length(breaksList)),
                     show_rownames = FALSE,
                     ... = ...)
  return(plot.pheatmap)
}


#' heatmap_differential
#'
#' plot differentially accessible peaks for a given comparison
#'
#' @param results cinaR result object
#' @param comparison these are created by cinaR from `contrasts` user provided. If not selected the first comparison will be shown!
#' @param ... additional arguments for heatmap function, for more info `?pheatmap`
#' @return ggplot object
#' @examples
#' \donttest{
#' library(cinaR)
#' data(atac_seq_consensus_bm) # calls 'bed'
#'
#' # a vector for comparing the examples
#' contrasts <- sapply(strsplit(colnames(bed), split = "-", fixed = TRUE),
#'                     function(x){x[1]})[4:25]
#'
#' results <- cinaR(bed, contrasts, reference.genome = "mm10")
#'
#' heatmap_differential(results)
#' }
#' @export
heatmap_differential <- function(results, comparison = NULL, ...){

  available.comparisons <- show_comparisons(results)

  if(is.null(comparison)){
    comparison <- available.comparisons[1]
    warning(paste0("'comparison' is not set so '", comparison, "' will be used!"))
  } else {
    # check if the comparison is available
    if(!comparison %in% available.comparisons){
      stop("'comparison' is not available. Please select one of these:\n",
           paste0("[",c(1:length(available.comparisons)),"] ", available.comparisons,collapse = "\n"))
    }
  }

  # if enrichment not run!
  if(!is.null(results[["cp"]])){
    cp <- results[["cp"]]
    da.peaks <- results$DA.peaks[[comparison]]
  } else { # if run
    cp <- results[["DA.results"]][["cp"]]
    da.peaks <- results$DA.results$DA.peaks[[comparison]]
  }

  cp <- stats::na.omit(cp)

  # select the peaks from contrast
  cp <- cp[da.peaks$Row.names,]

  # normalize peak distributions
  # utils function
  mat.heatmap <- scale_rows(cp)

  breaksList <- seq(min(mat.heatmap), max(mat.heatmap), by = .01)
  plot.pheatmap <- pheatmap::pheatmap(mat.heatmap, scale = "none",
                                      color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n=7,"RdBu")))(length(breaksList)),
                                      show_rownames = FALSE,
                                      ... = ...)
  return(plot.pheatmap)
}


