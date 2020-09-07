#' cinaR
#'
#' Runs DA pipeline and makes it ready for enrichment analyses
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
#' @return DApeaks returns DA peaks
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
           reference.genome = NULL) {

    if(length(contrasts) != (ncol(consensus.peaks)-3)){
      stop("Length of 'contrasts' must be equal to number of samples in 'consensus.peaks'")
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


    # edgeR, limma-voom, DEseq 2
    if (DA.choice %in% c(1:4)) {
      DA.results <- differentialAnalyses(cp = final.peaks,
                                         contrasts = contrasts,
                                         DA.choice = DA.choice,
                                         DA.fdr.threshold = DA.fdr.threshold,
                                         DA.lfc.threshold = DA.lfc.threshold,
                                         save.DA.peaks = save.DA.peaks)
    } else {
      stop (
        "DA.choice must be one of the followings.\n(1) edgeR, (2) limma-voom, (3) limma-trend, (4) DEseq2"
      )
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
#'
#' @return Normalized consensus peaks
#'
#' @export
normalizeConsensus <-
  function(cp, norm.method = "cpm") {
    if (norm.method == "cpm") {
      cp.norm <- edgeR::cpm(cp, log = T)
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

    bed <- as.data.frame(do.call(rbind, strsplit(rownames(cp), "_", fixed = T)))
    colnames(bed) <- c("CHR", "Start", "End")
    bed.GRanges <- GenomicRanges::GRanges(bed)

    if (is.null(reference.genome)) {
      warning("'reference.genome' is not set, therefore hg38 will be used!")
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
      genome <- annotables::grch38
    } else if (reference.genome == "hg38") {
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
      genome <- annotables::grch38
    } else if (reference.genome == "hg19") {
      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
      genome <- annotables::grch37
    } else if (reference.genome == "mm10") {
      txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
      genome <- annotables::grcm38
    } else {
      stop ("reference.genome should be 'hg38', 'hg19' or 'mm10'")
    }

    # annotate peaks
    annoPeaks <- ChIPseeker::annotatePeak(bed.GRanges, TxDb = txdb)

    if (show.annotation.pie) {
      ChIPseeker::plotAnnoPie(annoPeaks)
    }

    annoPeaks.anno <- annoPeaks@anno
    entrezids <- unique(annoPeaks.anno$geneId)

    # entrez to gene name mapping
    entrez2gene <- base::subset(genome, genome$entrez %in% entrezids, select = c('entrez', 'symbol'))

    # Match to each annotation dataframe
    m <- match(annoPeaks.anno$geneId, entrez2gene$entrez)
    annoPeaks.anno$gene_name <- entrez2gene$symbol[m]

    return(cbind(annoPeaks.anno, cp))
  }

#' Differential Analyses
#'
#' Runs differential analyses pipeline of choice on consensus peaks
#'
#' @param cp Consensus peaks
#' @param DA.choice determines which pipeline to run:
#' (1) edgeR, (2) limma-voom, (3) limma-trend, (4) DEseq2
#' @param contrasts user-defined contrasts for comparing samples
#' @param DA.fdr.threshold fdr cut-off for differential analyses
#' @param DA.lfc.threshold log-fold change cutoff for differential analyses
#' @param save.DA.peaks saves differentially accessible peaks to an excel file
#' @param DA.peaks.path the path which the excel file of the DA peaks will be saved,
#' if not set it will be saved to current directory.
#' @return DApeaks returns DA peaks
#'
#' @export
differentialAnalyses <- function(cp, contrasts, DA.choice,
                                 DA.fdr.threshold, DA.lfc.threshold,
                                 save.DA.peaks, DA.peaks.path) {
  cp.meta <- cp[, 1:15]
  cp.metaless <- cp[, 16:ncol(cp)]

  if (DA.choice %in% c(1:3)) { ## edgeR, limma-voom, limma-trend

    y <- edgeR::DGEList(counts = cp.metaless, group = contrasts)

    # Calculate normalization factors for library sizes with TMM
    y <- edgeR::calcNormFactors(y, method = "TMM")

    # Add intercept term for multiple comparisons
    design <- stats::model.matrix(~0 + contrasts)
    rownames(design) <- colnames(cp.metaless)
    colnames(design) <- levels(y$samples$group)

    # Estimate dispersion for genes with Bayesian Shrinkage
    cat(">> Estimating dispersion...\n")
    y <- edgeR::estimateDisp(y, design)

    # Fit the model
    cat(">> Fitting GLM...\n")
    fit.glm <- edgeR::glmQLFit(y, design)


    # Create contrasts for all comparisons
    combs <- utils::combn(colnames(design), 2)

    contrast.names <- apply(combs, 2,function(x){paste(x, collapse = "_")})
    cc <- apply(combs, 2,
          function(x){
            paste0(paste(x, collapse = "_"), "=", x[1], "-", x[2])
          })

    # to avoid the message in R CMD check!
    ccc <- NULL

    # create contrasts for package
    eval(parse(text = paste0("ccc <- limma::makeContrasts(",
                             paste(cc, collapse = ","),
                             ",levels = fit.glm$design)")))

    # Create DE gene list for differentially accessible peaks
    DA.peaks <- list()

    if(DA.choice == 1){ ## edgeR
      cat(">> Method: edgeR\n\tFDR:", DA.fdr.threshold, "& abs(logFC)<", DA.lfc.threshold, "\n")
      for (i in seq_len(ncol(ccc))){
        contrast.name <- colnames(ccc)[i]
        qlf <- edgeR::glmQLFTest(fit.glm,contrast = ccc[,i])
        # plotMD(qlf, main = contrast.name, p.value = 0.1)
        top.table <- edgeR::topTags(qlf, n = Inf, p.value = DA.fdr.threshold)$table
        top.table <- merge(cp.meta, top.table, by = 0)
        DA.peaks[[contrast.name]] <- top.table[abs(top.table$logFC) >= DA.lfc.threshold,]
      }
    } else if (DA.choice == 2) { ## limma-voom
        cat(">> Method: limma-voom\n\tFDR:", DA.fdr.threshold, "& abs(logFC)<", DA.lfc.threshold, "\n")
        v <- limma::voom(cp.metaless, design, plot=F)
        fit.voom <- limma::lmFit(v, design)
        fit.voom2 <- limma::eBayes(limma::contrasts.fit(fit.voom, ccc))
        # summary(decideTests(fit.voom2, method="separate", lfc = 0, p.value = 0.1))

        for (i in seq_len(ncol(ccc))){
          contrast.name <- colnames(ccc)[i]
          top.table <- limma::topTable(fit.voom2, coef = colnames(ccc)[i],
                                       p.value = DA.fdr.threshold, lfc = DA.lfc.threshold,
                                       number = Inf)
          top.table <- merge(cp.meta, top.table, by = 0)
          DA.peaks[[contrast.name]] <- top.table
        }
    } else if (DA.choice == 3){ ## limma-trend
        cat(">> Method: limma-trend\n\tFDR:", DA.fdr.threshold, "& abs(logFC)<", DA.lfc.threshold, "\n")
        fit.trend <- limma::lmFit(cp.metaless, design)
        fit.trend2 <- limma::eBayes(limma::contrasts.fit(fit.trend, ccc),
                                    trend = T)
        # summary(decideTests(fit.trend2, method="separate", lfc = 0, p.value = 0.1))

        for (i in seq_len(ncol(ccc))){
          contrast.name <- colnames(ccc)[i]
          top.table <- limma::topTable(fit.trend2, coef = colnames(ccc)[i],
                                p.value = DA.fdr.threshold,
                                lfc = DA.lfc.threshold,
                                number = Inf)
          top.table <- merge(cp.meta, top.table, by = 0)
          DA.peaks[[contrast.name]] <- top.table
        }
    } else if (DA.choice == 4) { ## DEseq2
        cat(">> Method: DEseq2\n\tFDR:", DA.fdr.threshold, "& abs(logFC)<", DA.lfc.threshold, "\n")

        # Create your desired groups
        group <- contrasts

        # Assign each sample to its group
        colData <- as.data.frame(cbind(colnames(cp.metaless), group))
        colnames(colData)  = c("sample", "groups")

        # Create DEseq Object
        dds <- DESeq2::DESeqDataSetFromMatrix(countData = cp.metaless, colData = colData, design = ~ groups)

        dds = DESeq2::DESeq(dds, parallel = T)

        # Create DE gene list for DESeq2

        for (i in seq_len(ncol(ccc))){
          contrast.name <- colnames(ccc)[i]
          DEseq.contrast <- rownames(ccc)[ccc[,i] != 0]
          res <- DESeq2::results(dds, c("groups", DEseq.contrast[2], DEseq.contrast[1]), parallel = T, tidy = T)
          rownames(res) <- res$row
          res.ordered <- res[order(res$pvalue),]
          res.significant <- subset(res.ordered, padj <= DA.fdr.threshold)
          res.significant <- subset(res.significant, abs(log2FoldChange) >= DA.lfc.threshold)
          res.significant <- merge(cp.meta, res.significant, by = 0)

          DA.peaks[[contrast.name]] <- res.significant
        }
    }
  } # end-if

  cat(">> DA peaks are found!\n")

  if(save.DA.peaks){
    if (is.null(DA.peaks.path)){
      cat(">> Saving DA peaks to current directory as DApeaks.xlsx...\n")
      writexl::write_xlsx(x = DA.peaks, path = "./DApeaks.xlsx")
    } else {
      cat(paste0(">> Saving DA peaks to ", DA.peaks.path,"...\n"))
      writexl::write_xlsx(x = DA.peaks, path = DA.peaks.path)
    }
  }

  return(DA.peaks)
}
