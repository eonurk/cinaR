
#' HPEA
#'
#' Hyper-geometric p-value enrichment analyses, looking for over-representation of a set of genes on given pathways.
#'
#' @param genes DA gene names to be checked if they are over-represented or not.
#' @param background.genes.size number of background genes for hyper-geometric p-value
#' calculations. Default is 20,000.
#' @param geneset Pathways to be used in enrichment analyses. If not set vp2008 (Chaussabel, 2008)
#' immune modules will be used. This can be set to any geneset using `read.gmt` function from `qusage`
#' package. Different modules are available: https://www.gsea-msigdb.org/gsea/downloads.jsp.
#' @return data.frame, list of pathways and their enrichment (adjusted) p-values.
#' @examples
#' \donttest{
#' library(cinaR)
#'
#' data("VP2008")
#' genes.to.test <- vp2008[[1]][1:10]
#' HPEA(genes.to.test,vp2008, background.genes.size = 20e3)
#' }
#'
#' @export
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
             lower.tail = FALSE,
             log.p = FALSE)

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
#' Gene set enrichment analyses, runs 'fgsea' package implementation with preset values.
#'
#' @param genes DA gene names to be checked if they are over-represented or not.
#' @param geneset Pathways to be used in enrichment analyses. If not set vp2008 (Chaussabel, 2008)
#' immune modules will be used. This can be set to any geneset using `read.gmt` function from `qusage`
#' package. Different modules are available: https://www.gsea-msigdb.org/gsea/downloads.jsp.
#' @return data.frame, list of pathways and their enrichment (adjusted) p-values.
#'
#' @examples
#' \donttest{
#' library(cinaR)
#' library(fgsea)
#' data(examplePathways)
#' data(exampleRanks)
#' GSEA(exampleRanks, examplePathways)
#' }
#' @references
#' G. Korotkevich, V. Sukhov, A. Sergushichev. Fast gene set enrichment analysis. bioRxiv (2019),
#' doi:10.1101/060012
#' @export
GSEA <- function(genes, geneset) {
  fgseaRes <- fgsea:: fgsea(pathways = geneset,
                            stats    = genes,
                            eps      = 0.0,
                            minSize  = 15,
                            maxSize  = 500)

  # Unlist the genes inside fgsea object to make it saveable.
  fgseaRes$leadingEdge <- vapply(fgseaRes$leadingEdge,function(x){paste(x,collapse = ",")},
                                 FUN.VALUE = "c")

  return(fgseaRes)
}

#' run_enrichment
#'
#' This function is run, if the enrichment pipeline wants to be called afterwards.
#' Setting reference genome to the same genome which cinaR was run should be given to this function!
#'
#' @param results list, DA peaks list for different contrasts
#' @param enrichment.method There are two methodologies for enrichment analyses,
#' Hyper-geometric p-value (HPEA) or Geneset Enrichment Analyses (GSEA).
#' @param experiment.type The type of experiment either set to "ATAC-Seq" or "RNA-Seq"
#' @param enrichment.FDR.cutoff FDR cut-off for enriched terms, p-values
#' are corrected by Benjamini-Hochberg procedure
#' @param reference.genome genome of interested species. It should be 'hg38', 'hg19' or 'mm10'.
#' @param background.genes.size number of background genes for hyper-geometric p-value
#' calculations. Default is 20,000.
#' @param geneset Pathways to be used in enrichment analyses. If not set vp2008 (Chaussabel, 2008)
#' immune modules will be used. This can be set to any geneset using `read.gmt` function from `qusage`
#' package. Different modules are available: https://www.gsea-msigdb.org/gsea/downloads.jsp.
#' @param verbose prints messages through running the pipeline
#' @return list, enrichment analyses results along with corresponding differential analyses outcomes
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
#' results <- cinaR(bed, contrasts, reference.genome = "mm10", run.enrichment = FALSE)
#'
#' results_with_enrichment <- run_enrichment(results, reference.genome = "mm10")
#' }
#' @export
run_enrichment <- function (
  results,
  geneset = NULL,
  experiment.type = "ATAC-Seq",
  reference.genome = NULL,
  enrichment.method = NULL,
  enrichment.FDR.cutoff = 1,
  background.genes.size = 20e3,
  verbose = TRUE
) {


  # Printing function
  verbosePrint <- verboseFn(verbose)

  # adjustment due to change in pipeline
  results <- results[["DA.peaks"]]

  # If no geneset is specified use VP2008
  if (is.null(geneset)) {
    verbosePrint(">> No `geneset` is specified so immune modules (Chaussabel, 2008) will be used!")
    geneset <- cinaR::vp2008
  }

  if (is.null(enrichment.method)) {
    verbosePrint(">> enrichment.method` is not selected. Hyper-geometric p-value (HPEA) will be used!")
    enrichment.method <- "HPEA"
  }

  if (reference.genome == "mm10") {

    # mm10 map
    mice2humanMap <- cinaR::grcm38

    # If it's an RNA-seq experiment and ensembl ids are used!
    if(experiment.type == "RNA-Seq" & grepl("ENSMUS", results[[1]][1,1], fixed = TRUE) ){

      results <- lapply(results, function(x) {
        m <- match(x[, "gene_name"], mice2humanMap[, "ensgene"])
        mapped.genes <- mice2humanMap[, "HGNC.symbol"][m]

        x[, "gene_name"] <- mapped.genes
        verbosePrint(">> Mice ensembl ids are converted to human symbols!")
        return(x)
      })

    } else { # If ATAC-Seq or gene-symbol are given

      # Convert to human homologs for enrichment analyses
      results <- lapply(results, function(x) {

        m <- match(x[, "gene_name"], mice2humanMap[, "symbol"])
        mapped.genes <- mice2humanMap[, "HGNC.symbol"][m]

        ## Careful: some genes may not be mappable!
        ## x[,"gene_name"] <- ifelse(is.na(mapped.genes),
        ##                          x[,"gene_name"], mapped.genes)
        x[, "gene_name"] <- mapped.genes

        return(x)
      })

      verbosePrint(">> Mice gene symbols are converted to human symbols!")
    }


  } else if (reference.genome == "hg38" & experiment.type == "RNA-Seq"){

    ens2gene <- cinaR::grch38

    if(grepl("ENSG", results[[1]][1,1], fixed = TRUE)){
      results <- lapply(results, function(x) {
        m <- match(x[, "gene_name"], ens2gene[, "ensgene"])
        mapped.genes <- ens2gene[, "symbol"][m]

        x[, "gene_name"] <- mapped.genes
        if(nrow(x) == sum(is.na(x[,"gene_name"])) & nrow(x) > 5){
          warning("Please make sure you are using the correct `reference.genome`!")
        }
        return(x)
      })
      verbosePrint(">> Human ensembl ids are converted to symbols...")
    }
  } else if (reference.genome == "hg19" & experiment.type == "RNA-Seq"){
    ens2gene <- cinaR::grch37
    if(grepl("ENSG", results[[1]][1,1], fixed = TRUE)){
      results <- lapply(results, function(x) {
        m <- match(x[, "gene_name"], ens2gene[, "ensgene"])
        mapped.genes <- ens2gene[, "symbol"][m]

        x[, "gene_name"] <- mapped.genes
        return(x)
      })
      verbosePrint(">> Human ensembl ids are converted to symbols...")
    }
  }

  # remove NA rows
  results <- lapply(results, function(x){

      if (length(x) > 0){ # check empty DA peaks
        x <- x[!is.na(x[, "gene_name"]),]

        if (experiment.type == "ATAC-Seq"){
          x <- x[order(abs(x[,"distanceToTSS"])),]
          x <- x[!duplicated(x[,"gene_name"]),]
        }
      }
    }
  )

  if (enrichment.method == "HPEA") {
    enrichment.results <-
    lapply(results, function(x) {

      if(is.null(x)){
        enrichment.results <- list()
      } else{
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
      }
      return(enrichment.results)
    })


  } else if (enrichment.method == "GSEA") {
    enrichment.results <- lapply(results, function(x){

      if(length(x) > 0) { # check non DA peaks
        rank <- x[,"logFC"]
        names(rank) <- x[,"gene_name"]
        res <- as.data.frame(GSEA(genes = rank, geneset = geneset))
        res <- res[res[,"padj"] <= enrichment.FDR.cutoff,]

      } else{
        res <- data.frame()
      }

      return(res)
    })
  } else {
    stop("`enrichment.method` should be either `HPEA` or `GSEA`")
  } # end-if (enrichment.method)

  # remove empty enrichment results
  enrichment.results <- Filter(function(x){length(x) >0}, enrichment.results)

  if(length(enrichment.results) == 0){
    warning(">> There are no enriched pathways for any of the comparisons!
            Maybe less stringent `DA.fdr.threshold` or `DA.lfc.threshold` values?")
  }
  return(enrichment.results)
}
