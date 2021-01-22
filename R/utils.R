#' mouse2human
#'
#' Given the mice gene symbols, this function creates a map from mice to human
#' using biomaRt.
#'
#' @param genes mice genes to be converted to human
#' @return returns a mapping from mouse to human
#'
#'
#' @examples
#' \donttest{
#'
#' mouse.genes <- c("Gfap", "Gzmb", "Il1b")
#' map <- mouse2human(mouse.genes)
#'
#' }
#' @export
mouse2human <- function(genes) {
  message(">> Converting mouse genes to human...\n")
  human = biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = biomaRt::useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")

  genesV2 = biomaRt::getLDS(
    attributes = c("mgi_symbol"),
    filters = "mgi_symbol",
    values = genes,
    mart = mouse,
    attributesL = c("hgnc_symbol"),
    martL = human,
    uniqueRows = TRUE
  )

  message(">> Mouse to human mapping is created for all genes!\n")
  # Return map
  return(genesV2)
}

#' scale_rows
#'
#' Normalize (z-score) rows of a matrix
#'
#' @param x a matrix, possibly containing gene by samples
#' @return Row-normalized matrix
#'
#'
#' @examples
#' \donttest{
#' library(cinaR)
#' data(atac_seq_consensus_bm) # calls 'bed'
#' bed.row.normalized <- scale_rows(bed[,c(4:25)])
#' head(bed.row.normalized)
#' }
#' @export
scale_rows <- function (x){
  m <- apply(x, 1, mean, na.rm = TRUE)
  s <- apply(x, 1, stats::sd, na.rm = TRUE)
  return((x - m) / s)
}


#' show_comparisons
#'
#' returns the names of the created comparisons
#'
#' @param results output of the cinaR
#' @return comparisons created
#' @export
show_comparisons <- function (results) {
  # if enrichment not run!
  if(!is.null(results[["cp"]])){
    return(names(results$DA.peaks))
  } else { # if run
    return(names(results$DA.results$DA.peaks))
  }
}

#' verboseFn
#'
#' returns a printing function to be used with in the script
#' @param verbose boolean, determines whether the output going be printed or not
#' @return print function
#'
#' @export
verboseFn <- function(verbose){
  if(verbose){
    return(function(...) {
      message(...)
    })
  } else {
    return(function(...){
      return(invisible(NULL))
    })
  }
}


