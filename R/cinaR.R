#' cinaR
#'
#' Runs differential analyses and enrichment pipelines
#'
#' @param consensus.peaks bed formatted consensus peak matrix: CHR, START, STOP and raw peak counts (peaks by 3+samples)
#' @param contrasts user-defined contrasts for comparing samples
#' @param DA.choice determines which pipeline to run:
#' (1) edgeR, (2) limma-voom, (3) limma-trend, (4) DEseq2.
#' Note: Use limma-trend if consensus peaks are already normalized, otherwise use other methods.
#' @param norm.method normalization method for consensus peaks
#' @param filter.method filtering method for low expressed peaks
#' @param TSS.threshold Distance to transcription start site in base-pairs. Default set to 50,000.
#' @param library.threshold number of libraries a peak occurs so that it is not filtered default set to 2
#' @param cpm.threshold count per million threshold for not to filter a peak
#' @param show.annotation.pie shows the annotation pie chart produced with ChipSeeker
#' @param reference.genome genome of interested species. It should be 'hg38', 'hg19' or 'mm10'.
#' @param DA.fdr.threshold fdr cut-off for differential analyses
#' @param DA.lfc.threshold log-fold change cutoff for differential analyses
#' @param save.DA.peaks saves differentially accessible peaks to an excel file
#' @param DA.peaks.path the path which the excel file of the DA peaks will be saved,
#' if not set it will be saved to current directory.
#' @param batch.correction logical, if set will run unsupervised batch correction
#' via sva (default) or if the batch information is known `batch.information`
#' argument should be provided by user.
#' @param batch.information character vector, given by user.
#' @param run.enrichment logical, turns off enrichment pipeline
#' @param enrichment.method There are two methodologies for enrichment analyses,
#' Hyper-geometric p-value (HPEA) or Geneset Enrichment Analyses (GSEA).
#' @param enrichment.FDR.cutoff FDR cut-off for enriched terms, p-values
#' are corrected by Benjamini-Hochberg procedure
#' @param background.genes.size number of background genes for hyper-geometric p-value
#' calculations. Default is 20,000.
#' @param geneset Pathways to be used in enrichment analyses. If not set vp2008 (Chaussabel, 2008)
#' immune modules will be used. This can be set to any geneset using `read.gmt` function from `qusage`
#' package. Different modules are available: https://www.gsea-msigdb.org/gsea/downloads.jsp.
#'
#' @examples
#' \dontrun{
#' data(atac_seq_consensus_bm) # calls 'bed'
#'
#' # a vector for comparing the examples
#' contrasts <- sapply(strsplit(colnames(bed), split = "-", fixed = TRUE),
#'                     function(x){x[1]})[4:25]
#'
#' DA.results <- cinaR(bed, contrasts, reference.genome = "mm10")
#' }
#'
#'
#' @return returns differentially accessible peaks
#'
#' @export
cinaR <-
  function(consensus.peaks,
           contrasts,
           DA.choice = 1,
           DA.fdr.threshold = 0.05,
           DA.lfc.threshold = 0,
           save.DA.peaks = F,
           DA.peaks.path = NULL,
           norm.method = "cpm",
           filter.method = "custom",
           library.threshold = 2,
           cpm.threshold = 1,
           TSS.threshold = 50e3,
           show.annotation.pie = F,
           reference.genome = NULL,
           batch.correction = FALSE,
           batch.information = NULL,
           run.enrichment = TRUE,
           enrichment.method = NULL,
           enrichment.FDR.cutoff = 1,
           background.genes.size = 20e3,
           geneset = NULL) {

    if (length(contrasts) != (ncol(consensus.peaks) - 3)) {
      stop("Length of 'contrasts' must be equal to number of samples in 'consensus.peaks'")
    }

    if (is.null(reference.genome)) {
      warning("'reference.genome' is not set, therefore hg38 will be used!")
      reference.genome <- "hg38"
    }

    # collapse chr, start, end and make them rownames
    cp.rownames <-
      apply(consensus.peaks[, 1:3], 1, function(x) {
        paste0(trimws(x), collapse = "_")
      })
    cp <- consensus.peaks[, -c(1:3)]
    rownames(cp) <- cp.rownames

    # filter low expressed peaks
    cp.filtered <-
      filterConsensus(cp, library.threshold = library.threshold, cpm.threshold = cpm.threshold)


    # annotate the peaks to the closest TSS
    cp.filtered.annotated <- annotatePeaks(cp.filtered,
                                           reference.genome = reference.genome,
                                           show.annotation.pie = show.annotation.pie)

    # filter distance to TSS
    final.peaks <-
      cp.filtered.annotated[abs(cp.filtered.annotated$distanceToTSS) <= TSS.threshold, ]


    if (!is.null(enrichment.method)) {
      if (enrichment.method == "GSEA" & run.enrichment == TRUE) {
        message(
          ">> Setting `DA.fdr.threshold = 1` and `DA.lfc.threshold = 0`
              since GSEA is chosen for enrichment method!"
        )

        DA.fdr.threshold <- 1
        DA.lfc.threshold <- 0
      }
    }

    # edgeR, limma-voom, DEseq 2
    if (DA.choice %in% c(1:4)) {
      DA.results <- differentialAnalyses(
        final.peaks = final.peaks,
        contrasts = contrasts,
        DA.choice = DA.choice,
        DA.fdr.threshold = DA.fdr.threshold,
        DA.lfc.threshold = DA.lfc.threshold,
        save.DA.peaks = save.DA.peaks,
        batch.correction = batch.correction,
        batch.information = batch.information
      )
    } else {
      stop (
        "DA.choice must be one of the followings.
        (1) edgeR, (2) limma-voom, (3) limma-trend, (4) DEseq2"
      )
    }

    if (run.enrichment) {
      enrichment.results <-
        run_enrichment(
          results = DA.results,
          geneset = geneset,
          reference.genome = reference.genome,
          enrichment.method = enrichment.method,
          enrichment.FDR.cutoff = enrichment.FDR.cutoff,
          background.genes.size = background.genes.size
        )
      cat(">> Enrichment results are ready...\n")
      cat(">> Done!\n")
      return(list(DA.results = DA.results,
                  Enrichment.Results = enrichment.results))
    }
    return(DA.results)
  }

#' filterConsensus
#'
#' Filters lowly expressed peaks from down-stream analyses
#'
#' @importFrom edgeR cpm filterByExpr
#' @param cp bed formatted consensus peak matrix: CHR, START, STOP and raw peak counts (peaks by 3+samples)
#' @param filter.method filtering method for low expressed peaks
#' @param library.threshold number of libraries a peak occurs so that it is not filtered default set to 2
#' @param cpm.threshold count per million threshold for not to filter a peak
#'
#' @return DApeaks returns DA peaks
#'
#' @export
filterConsensus <-
  function(cp,
           filter.method = "custom",
           library.threshold = 2,
           cpm.threshold = 1) {
    if (filter.method == "custom") {
      cp.filtered <-
        cp[rowSums(edgeR::cpm(cp) >= cpm.threshold) >= library.threshold, ]
    } else if (filter.method == "edgeR") {
      cp.filtered <- edgeR::filterByExpr(cp)
    } else {
      stop("filter.method should be either 'custom' or 'edgeR'")
    }
  }


#' normalizeConsensus
#'
#' Normalizes consensus peak using different methods
#'
#' @param cp bed formatted consensus peak matrix: CHR, START, STOP and raw peak counts (peaks by 3+samples)
#' @param norm.method normalization method for consensus peaks
#' @param log.option logical, log option for cpm function in edgeR
#' @return Normalized consensus peaks
#'
#' @export
normalizeConsensus <-
  function(cp, norm.method = "cpm", log.option = FALSE) {
    if (norm.method == "cpm") {
      cp.norm <- edgeR::cpm(cp, log = log.option, prior.count = 5)
    } else if (norm.method == "quantile") {
      cp.norm <- preprocessCore::normalize.quantiles(cp)
    } else {
      stop("Wrong normalization method, it must be either 'cpm' or 'quantile'")
    }
    return(cp.norm)
  }

#' annotatePeaks
#'
#' Runs DA pipeline and makes it ready for enrichment analyses
#'
#' @param cp bed formatted consensus peak matrix: CHR, START, STOP and raw peak counts (peaks by 3+samples)
#' @param reference.genome genome of interested species. It should be 'hg38', 'hg19' or 'mm10'.
#' @param show.annotation.pie shows the annotation pie chart produced with ChipSeeker
#'
#' @return DApeaks returns DA peaks
#'
#' @export
annotatePeaks <-
  function(cp,
           reference.genome,
           show.annotation.pie = F) {
    bed <-
      as.data.frame(do.call(rbind, strsplit(rownames(cp), "_", fixed = T)))
    colnames(bed) <- c("CHR", "Start", "End")
    bed.GRanges <- GenomicRanges::GRanges(bed)


    if (reference.genome == "hg38") {
      if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
        stop(
          "Package \"TxDb.Hsapiens.UCSC.hg38.knownGene\" needed for this
             function to work. Please install it.",
          call. = FALSE
        )
      }
      txdb <-
        TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
      genome <- cinaR::grch38
      reference.genome <- "hg38"
    } else if (reference.genome == "hg19") {
      if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE)) {
        stop(
          "Package \"TxDb.Hsapiens.UCSC.hg19.knownGene\" needed for this
             function to work. Please install it.",
          call. = FALSE
        )
      }
      txdb <-
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
      genome <- cinaR::grch37
    } else if (reference.genome == "mm10") {
      if (!requireNamespace("TxDb.Mmusculus.UCSC.mm10.knownGene", quietly = TRUE)) {
        stop(
          "Package \"TxDb.Mmusculus.UCSC.mm10.knownGene\" needed for this
             function to work. Please install it.",
          call. = FALSE
        )
      }
      txdb <-
        TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
      genome <- cinaR::grcm38
    } else {
      stop ("reference.genome should be 'hg38', 'hg19' or 'mm10'")
    }

    # annotate peaks
    annoPeaks <- ChIPseeker::annotatePeak(bed.GRanges, TxDb = txdb, assignGenomicAnnotation = T)

    if (show.annotation.pie) {
      ChIPseeker::plotAnnoPie(annoPeaks)
    }

    annoPeaks.anno <- annoPeaks@anno
    entrezids <- unique(annoPeaks.anno$geneId)

    # entrez to gene name mapping
    entrez2gene <-
      base::subset(genome,
                   genome$entrez %in% entrezids,
                   select = c('entrez', 'symbol'))

    # Match to each annotation dataframe
    m <- match(annoPeaks.anno$geneId, entrez2gene$entrez)
    annoPeaks.anno$gene_name <- entrez2gene$symbol[m]

    return(cbind(annoPeaks.anno, cp))
  }

#' Differential Analyses
#'
#' Runs differential analyses pipeline of choice on consensus peaks
#'
#' @param final.peaks Annotated Consensus peaks
#' @param DA.choice determines which pipeline to run:
#' (1) edgeR, (2) limma-voom, (3) limma-trend, (4) DEseq2
#' @param contrasts user-defined contrasts for comparing samples
#' @param DA.fdr.threshold fdr cut-off for differential analyses
#' @param DA.lfc.threshold log-fold change cutoff for differential analyses
#' @param save.DA.peaks logical, saves differentially accessible peaks to an excel file
#' @param DA.peaks.path the path which the excel file of the DA peaks will be saved,
#' if not set it will be saved to current directory.
#' @param batch.correction logical, if set will run unsupervised batch correction
#' via sva (default) or if the batch information is known `batch.information`
#' argument should be provided by user.
#' @param batch.information character vector, given by user.
#'
#' @return returns consensus peaks (batch corrected version if enabled) and DA peaks
#'
#' @export
differentialAnalyses <- function(final.peaks,
                                 contrasts,
                                 DA.choice,
                                 DA.fdr.threshold,
                                 DA.lfc.threshold,
                                 save.DA.peaks,
                                 DA.peaks.path,
                                 batch.correction,
                                 batch.information) {

  # silence build notes
  log2FoldChange <- padj <- NULL

  cp.meta <- final.peaks[, 1:15]
  cp.metaless <- final.peaks[, 16:ncol(final.peaks)]

  design <- stats::model.matrix(~ 0 + contrasts)

  if (batch.correction) {
    if (is.null(batch.information)) {
      ## First normalize the consensus peaks to avoid detecting the effects
      ## confounding from library size as Michael Love and Jeff Leek suggests
      ## in this thread:
      cat(">> Running SVA for batch correction...\n")

      cp.metaless.normalized <- normalizeConsensus(cp.metaless)
      mod  <- stats::model.matrix(~ 0 + contrasts)
      mod0 <- cbind(rep(1, length(contrasts)))

      # calculate the batch effects
      sva.res <-
        sva::svaseq(cp.metaless.normalized, mod, mod0)

      # batch effect additional matrix
      add.batch <- sva.res$sv

      # make the colnames prettier just for fun
      colnames(add.batch) <- paste0("SV", c(1:ncol(add.batch)))

      # add it to the design matrix
      design <-
        cbind(design, add.batch)

      # batch corrected consensus peaks created for PCA/Heatmaps
      cp.batch.corrected <- limma::removeBatchEffect(cp.metaless, covariates = add.batch)

    } else {
      # if there is batch information available
      cat(">> Adding batch information to design matrix...\n")

      if (nrow(design) != length(batch.information)) {
        stop("Number of samples and `batch.information` should be same length!")
      }

      design <- cbind(design, BatchInfo = batch.information)

      # batch corrected consensus peaks created for PCA/Heatmaps
      cp.batch.corrected <- limma::removeBatchEffect(cp.metaless, batch = batch.information)
    }
  }

  # Add intercept term for multiple comparisons

  rownames(design) <- colnames(cp.metaless)
  colnames(design) <- gsub("contrasts", "", colnames(design))

  # Create contrasts for all comparisons
  combs <-
    utils::combn(colnames(design)[1:length(unique(contrasts))], 2)

  contrast.names <-
    apply(combs, 2, function(x) {
      paste(x, collapse = "_")
    })
  cc <- apply(combs, 2,
              function(x) {
                paste0(paste(x, collapse = "_"), "=", x[1], "-", x[2])
              })

  # to avoid the message in R CMD check!
  ccc <- NULL

  # create contrasts to be compared
  eval(parse(
    text = paste0(
      "ccc <- limma::makeContrasts(",
      paste(cc, collapse = ","),
      ",levels = colnames(design))"
    )
  ))

  # Create DE gene list for differentially accessible peaks
  DA.peaks <- list()

  if (DA.choice == 1) {
    ## edgeR

    cat(
      ">> Method: edgeR\n\tFDR:",
      DA.fdr.threshold,
      "& abs(logFC)<",
      DA.lfc.threshold,
      "\n"
    )

    y <- edgeR::DGEList(counts = cp.metaless, group = contrasts)

    # Calculate normalization factors for library sizes with TMM
    y <- edgeR::calcNormFactors(y, method = "TMM")

    # Estimate dispersion for genes with Bayesian Shrinkage
    cat(">> Estimating dispersion...\n")
    y <- edgeR::estimateDisp(y, design)

    # Fit the model
    cat(">> Fitting GLM...\n")
    fit.glm <- edgeR::glmQLFit(y, design)

    for (i in seq_len(ncol(ccc))) {
      contrast.name <- colnames(ccc)[i]
      qlf <- edgeR::glmQLFTest(fit.glm, contrast = ccc[, i])
      # plotMD(qlf, main = contrast.name, p.value = 0.1)
      top.table <-
        edgeR::topTags(qlf, n = Inf, p.value = DA.fdr.threshold)$table
      top.table <- merge(cp.meta, top.table, by = 0)

      # ifelse does not return the dataframe for some reason,
      # therefore, implemented this check explicitly
      if (nrow(top.table) > 0) {
        top.table <- top.table[abs(top.table$logFC) >= DA.lfc.threshold,]
        # Refactor to uniformize DA results
        top.table <- top.table[, c(1:17, 21)]
        DA.peaks[[contrast.name]] <- top.table
      } else {
        DA.peaks[[contrast.name]] <- list()
      }
    }
  } else if (DA.choice == 2) {
    ## limma-voom
    cat(
      ">> Method: limma-voom\n\tFDR:",
      DA.fdr.threshold,
      "& abs(logFC)<",
      DA.lfc.threshold,
      "\n"
    )
    v <- limma::voom(cp.metaless, design, plot = F)
    fit.voom <- limma::lmFit(v, design)
    fit.voom2 <-
      limma::eBayes(limma::contrasts.fit(fit.voom, ccc))
    # summary(decideTests(fit.voom2, method="separate", lfc = 0, p.value = 0.1))

    for (i in seq_len(ncol(ccc))) {
      contrast.name <- colnames(ccc)[i]
      top.table <-
        limma::topTable(
          fit.voom2,
          coef = colnames(ccc)[i],
          p.value = DA.fdr.threshold,
          lfc = DA.lfc.threshold,
          number = Inf
        )
      top.table <- merge(cp.meta, top.table, by = 0)

      # Refactor to uniformize DA results
      top.table <- top.table[, c(1:17, 21)]
      colnames(top.table)[18] <- "FDR"

      if (nrow(top.table) > 0) {
        DA.peaks[[contrast.name]] <- top.table
      } else {
        DA.peaks[[contrast.name]] <- list()
      }

    }
  } else if (DA.choice == 3) {
    ## limma-trend
    cat(
      ">> Method: limma-trend\n\tFDR:",
      DA.fdr.threshold,
      "& abs(logFC)<",
      DA.lfc.threshold,
      "\n"
    )
    fit.trend <- limma::lmFit(cp.metaless, design)
    fit.trend2 <-
      limma::eBayes(limma::contrasts.fit(fit.trend, ccc),
                    trend = T)
    # summary(decideTests(fit.trend2, method="separate", lfc = 0, p.value = 0.1))

    for (i in seq_len(ncol(ccc))) {
      contrast.name <- colnames(ccc)[i]
      top.table <-
        limma::topTable(
          fit.trend2,
          coef = colnames(ccc)[i],
          p.value = DA.fdr.threshold,
          lfc = DA.lfc.threshold,
          number = Inf
        )
      top.table <- merge(cp.meta, top.table, by = 0)

      # Refactor to uniformize DA results
      top.table <- top.table[, c(1:17, 21)]
      colnames(top.table)[18] <- "FDR"

      if (nrow(top.table) > 0) {
        DA.peaks[[contrast.name]] <- top.table
      } else {
        DA.peaks[[contrast.name]] <- list()
      }
    }
  } else if (DA.choice == 4) {
    ## DEseq2
    cat(
      ">> Method: DEseq2\n\tFDR:",
      DA.fdr.threshold,
      "& abs(logFC)<",
      DA.lfc.threshold,
      "\n"
    )

    # Create your desired groups
    group <- contrasts

    # Assign each sample to its group
    colData <- as.data.frame(cbind(colnames(cp.metaless), group))
    colnames(colData)  = c("sample", "groups")

    # Create DEseq Object
    dds <-
      DESeq2::DESeqDataSetFromMatrix(countData = cp.metaless,
                                     colData = colData,
                                     design = ~ groups)

    dds = DESeq2::DESeq(dds, parallel = T)

    # Create DE gene list for DESeq2

    for (i in seq_len(ncol(ccc))) {
      contrast.name <- colnames(ccc)[i]
      DEseq.contrast <- rownames(ccc)[ccc[, i] != 0]
      res <-
        DESeq2::results(
          dds,
          c("groups", DEseq.contrast[2], DEseq.contrast[1]),
          parallel = T,
          tidy = T
        )
      rownames(res) <- res$row
      res.ordered <- res[order(res$pvalue), ]
      res.significant <-
        subset(res.ordered,
               padj <= DA.fdr.threshold &
                 abs(log2FoldChange) >= DA.lfc.threshold)
      res.significant <- merge(cp.meta, res.significant, by = 0)

      top.table <- res.significant[, c(1:16, 19, 23)]
      colnames(top.table)[c(17, 18)] <- c("logFC", "FDR")
      if (nrow(top.table) > 0) {
        DA.peaks[[contrast.name]] <- top.table
      } else {
        DA.peaks[[contrast.name]] <- list()
      }
    }
  }

  cat(">> DA peaks are found!\n")

  if (save.DA.peaks) {
    if (is.null(DA.peaks.path)) {
      cat(">> Saving DA peaks to current directory as DApeaks.xlsx...\n")
      writexl::write_xlsx(x = DA.peaks, path = "./DApeaks.xlsx")
    } else {
      cat(paste0(">> Saving DA peaks to ", DA.peaks.path, "...\n"))
      writexl::write_xlsx(x = DA.peaks, path = DA.peaks.path)
    }
  }

  if (batch.correction){
    return(list (cp = cp.batch.corrected, DA.peaks = DA.peaks))
  }
  return(list (cp = cp.metaless, DA.peaks = DA.peaks))
}
