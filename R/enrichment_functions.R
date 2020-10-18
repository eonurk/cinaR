
#' HPEA
#'
#' @param genes DA gene names to be checked if they are over-represented or not.
#' @param background.genes.size number of background genes for hyper-geometric p-value
#' calculations. Default is 20,000.
#' @param geneset Pathways to be used in enrichment analyses. If not set vp2008 (Chaussabel, 2008)
#' immune modules will be used. This can be set to any geneset using `read.gmt` function from `qusage`
#' package. Different modules are available: https://www.gsea-msigdb.org/gsea/downloads.jsp.
HPEA <- function(genes, geneset, background.genes.size) {
  res.hpea <- lapply(geneset, function(module) {
    # Module gene count
    k <- length(unique(module))

    # gene count
    n <- length(genes)

    # overlapped up-regulated genes with the module
    q <- sum(genes %in% module)


    # here we calculate the probability of having a bigger intersection
    # than the count of overlapping genes given the module size and the total gene count.
    # we substract 1 for removing the equality when the lower.tail = F, which changes P(X<x) to 1-P(X>=x).
    p.hyper <-
      stats::phyper(q - 1,
             k,
             background.genes.size - k,
             n,
             lower.tail = F,
             log.p = F)

    # take the overlapping genes for the modules
    og <- genes[genes %in% module]

    module.overlapping.ratio <- length(og) / k

    og <- ifelse(length(og) <= 0, "", paste(og, collapse = ","))


    return(
      data.frame(
        p.val = p.hyper,
        overlapping.genes = og,
        module.overlapping.ratio = module.overlapping.ratio
      )
    )
  })

  res.hpea <- do.call(rbind, res.hpea)

  res.hpea<- cbind(module.name = names(geneset), res.hpea)
  res.hpea[, "adj.p"] <- stats::p.adjust(res.hpea[, "p.val"], method = "BH")

  # sort according to adjusted p-values and then to p-values
  res.hpea <-
    res.hpea  [order(res.hpea[, "adj.p"], res.hpea[, "p.val"]), ]

  return(res.hpea)
}


#' GSEA
#'
#' @param genes DA gene names to be checked if they are over-represented or not.
#' @param geneset Pathways to be used in enrichment analyses. If not set vp2008 (Chaussabel, 2008)
#' immune modules will be used. This can be set to any geneset using `read.gmt` function from `qusage`
#' package. Different modules are available: https://www.gsea-msigdb.org/gsea/downloads.jsp.
#' TODO: add ... to parameters
#' @export
GSEA <- function(genes, geneset) {
  fgseaRes <- fgsea::fgsea(geneset, genes)

  # Unlist the genes inside fgsea object to make it saveable.
  fgseaRes$leadingEdge <- vapply(fgseaRes$leadingEdge,function(x){paste(x,collapse = ",")},
                                 FUN.VALUE = "c")

  return(fgseaRes)
}

#' run_enrichment
#'
#' @param results list, DA peaks list for different contrasts
#' @param enrichment.method There are two methodologies for enrichment analyses,
#' Hyper-geometric p-value (HPEA) or Geneset Enrichment Analyses (GSEA).
#' @param enrichment.FDR.cutoff FDR cut-off for enriched terms, p-values
#' are corrected by Benjamini-Hochberg procedure
#' @param reference.genome genome of interested species. It should be 'hg38', 'hg19' or 'mm10'.
#' @param background.genes.size number of background genes for hyper-geometric p-value
#' calculations. Default is 20,000.
#' @param geneset Pathways to be used in enrichment analyses. If not set vp2008 (Chaussabel, 2008)
#' immune modules will be used. This can be set to any geneset using `read.gmt` function from `qusage`
#' package. Different modules are available: https://www.gsea-msigdb.org/gsea/downloads.jsp.
run_enrichment <- function (
  results,
  geneset = NULL,
  reference.genome = NULL,
  enrichment.method = NULL,
  enrichment.FDR.cutoff = 0.05,
  background.genes.size = 20e3
) {
  # If no geneset is specified use VP2008
  if (is.null(geneset)) {
    message(">> No `geneset` is specified so immune modules (Chaussabel, 2008) will be used!")
    geneset <- cinaR::vp2008
  }


  if (is.null(reference.genome)){
    message(">> `reference.genome` is not specified, so DA genes are assumed to be human")
  }
  if (is.null(enrichment.method)) {
    message(">> enrichment.method` is not selected. Hyper-geometric p-value (HPEA) will be used!")
    enrichment.method <- "HPEA"
  }
  # Convert mice gene symbols to human's
  if (reference.genome == "mm10") {
    all.genes <- sapply(results, function(x) {
      x["gene_name"]
    })
    all.genes <- unlist(all.genes)
    all.genes <- unique(all.genes)

    # Create mapping for all mice genes
    mice2humanMap <- cinaR::mouse2human(genes = all.genes)

    # Convert to human homologs for enrichment analyses
    results <- lapply(results, function(x) {
      m <- match(x[, "gene_name"], mice2humanMap[, "MGI.symbol"])
      mapped.genes <- mice2humanMap[, "HGNC.symbol"][m]

      ## Careful: some genes may not be mappable!
      ## x[,"gene_name"] <- ifelse(is.na(mapped.genes),
      ##                          x[,"gene_name"], mapped.genes)
      x[, "gene_name"] <- mapped.genes
      cat(">> Human gene symbols are converted to mice!\n")
      return(x)
    })
  }

  results <- lapply(results, function(x){
    # remove NA rows
    x <- x[!is.na(x[, "gene_name"]),]
  })

  if (enrichment.method == "HPEA") {
    enrichment.results <-
    lapply(results, function(x) {
      opening.locs <- x[, "logFC"] > 0
      genes.opening <- x[ opening.locs, ]
      genes.closing <- x[!opening.locs, ]

      enrichment.opening <-
        HPEA(
          genes = genes.opening[, "gene_name"],
          geneset = geneset,
          background.genes.size = background.genes.size
        )

      enrichment.closing <-
        HPEA(
          genes = genes.closing[, "gene_name"],
          geneset = geneset,
          background.genes.size = background.genes.size
        )

      enrichment.results <-
        rbind(
          cbind(enrichment.opening, status = "Opening"),
          cbind(enrichment.closing, status = "Closing")
        )

      enrichment.results <-
        enrichment.results[enrichment.results[,"adj.p"] <= enrichment.FDR.cutoff,]

      rownames(enrichment.results) <- NULL

      return(enrichment.results)
    })
  } else if (enrichment.method == "GSEA") {
    enrichment.results <- lapply(results, function(x){
      rank <- x[,"logFC"]
      names(rank) <- x[,"gene_name"]
      res <- as.data.frame(GSEA(genes = rank, geneset = geneset))
      res <- res[res[,"padj"] <= enrichment.FDR.cutoff,]
      return(res)
    })
  } else {
    stop("`enrichment.method` should be either `HPEA` or `GSEA`")
  } # end-if (enrichment.method)

  return(enrichment.results)
}
