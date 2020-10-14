
#' dot_plot
#'
#' @param object cinaR result object
#' @param fdr.cutoff Pathways with smaller fdr values than the cut-off
#' will be shown as dots.
#' @param filter.pathways logical, it will filter the pathways from dot plot
#' with fdr values less than `fdr.cutoff`.
#' @return ggplot object
#'
#' @export
dot_plot <- function(object, fdr.cutoff = 0.1, filter.pathways = FALSE){

  source("R/color_values.R")

  if(is.null (object[["Enrichment.Results"]])){
    stop("Did you run the enrichment pipeline in cinaR? For more info ?cinaR")
  }

  object <- object[["Enrichment.Results"]]

  # add list names as a column
  df.plot <- dplyr::bind_rows(object, .id = "column_label")

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
