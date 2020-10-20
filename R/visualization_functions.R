
#' dot_plot
#'
#' @param results cinaR result object
#' @param fdr.cutoff Pathways with smaller fdr values than the cut-off
#' will be shown as dots.
#' @param filter.pathways logical, it will filter the pathways from dot plot
#' with fdr values less than `fdr.cutoff`.
#' @return ggplot object
#'
#' @export
dot_plot <- function(results, fdr.cutoff = 0.1, filter.pathways = FALSE){

  if(is.null (results[["Enrichment.Results"]])){
    stop("Did you run the enrichment pipeline in cinaR? For more info ?cinaR")
  }

  #silence build NOTES
  adj.p <- column_label <-
    log2FoldChange <- module.name <- padj <-
    prcomp <- status <- NULL

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
#' @return ggplot object, pca plot
#'
#' @export
pca_plot <- function(results, overlaid.info, sample.names = NULL, show.names = TRUE){

  #silence CRAN build NOTES
  PC1 <- PC2 <- NULL

  # extract consensus peaks
  cp <- results[["DA.results"]][["cp"]]

  if(is.null(sample.names)){
    sample.names <- colnames(cp)
  } else{
    if(length(sample.names) != ncol(cp)){
      stop("The length of `sample.names` should be equal to number of samples.")
    }
  }

  # consensus peaks are normalized before pca.
  # log option is set TRUE to have a better variance-stabilization.
  cp <- normalizeConsensus(cp, log.option = TRUE)


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

#' heatmap_plot
#'
#' plot most variable k peaks (default k = 100) among all samples
#'
#' @param results cinaR result object
#' @param heatmap.peak.count number of peaks to be plotted.
#' If number of peaks are less than k then all peaks will be used.
#' @param ... additional arguments for heatmap function, for more info `?pheatmap`
#' @return ggplot object, pca plot
#'
#' @export
heatmap_plot <- function(results, heatmap.peak.count = 100, ...){
  cp <- results[["DA.results"]][["cp"]]

  # consensus peaks are normalized before pca.
  # log option is set TRUE to have a better variance-stabilization.
  # NOTE: This is done because of the library depth differences, it does not
  # account for peak-peak differences
  cp <- normalizeConsensus(cp, log.option = TRUE)

  cp <- stats::na.omit(cp)

  # Remove possible na's from data and
  # reorder peak according to their standard deviation in decreasing order
  cp <- cp [rev(order(apply(cp, 1, stats::sd))),]

  mat.heatmap <- cp[1:min(nrow(cp),heatmap.peak.count),]

  # pheatmap normalization function
  scale_rows <- function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, stats::sd, na.rm = T)
    return((x - m) / s)
  }

  # normalize peak distributions
  mat.heatmap <- scale_rows(mat.heatmap)

  breaksList = seq(min(mat.heatmap), max(mat.heatmap), by = .01)
  plot.pheatmap <- pheatmap::pheatmap(mat.heatmap,
                     color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n=7,"RdBu")))(length(breaksList)),
                     show_rownames = F,
                     ... = ...)
  return(plot.pheatmap)
}


