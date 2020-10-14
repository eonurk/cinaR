#' mouse2human
#'
#' @param genes mice genes to be converted to human
#' @return returns a mapping from mouse to human
#'
#' @export
mouse2human <- function(genes) {
  cat(">> Converting mouse genes to human...\n")
  human = biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = biomaRt::useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")

  genesV2 = biomaRt::getLDS(
    attributes = c("mgi_symbol"),
    filters = "mgi_symbol",
    values = genes,
    mart = mouse,
    attributesL = c("hgnc_symbol"),
    martL = human,
    uniqueRows = T
  )

  cat(">> Mouse to human mapping is created for all genes!\n")
  # Return map
  return(genesV2)
}
