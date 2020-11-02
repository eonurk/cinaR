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
