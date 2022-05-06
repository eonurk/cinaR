#' cinaR
#'
#' Runs differential analyses and enrichment pipelines
#'
#' @param matrix either bed formatted consensus peak matrix (peaks by 3+samples)
#' CHR, START, STOP and raw peak counts OR count matrix (genes by 1+samples).
#' @param contrasts user-defined contrasts for comparing samples
#' @param experiment.type The type of experiment either set to "ATAC-Seq" or "RNA-Seq"
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
#' @param comparison.scheme either one-vs-one (OVO) or one-vs-all (OVA) comparisons.
#' @param save.DA.peaks saves differentially accessible peaks to an excel file
#' @param DA.peaks.path the path which the excel file of the DA peaks will be saved,
#' if not set it will be saved to current directory.
#' @param batch.correction logical, if set will run unsupervised batch correction
#' via sva (default) or if the batch information is known `batch.information`
#' argument should be provided by user.
#' @param batch.information character vector, given by user.
#' @param additional.covariates vector or data.frame, this parameter will be directly added to design
#' matrix before running the differential analyses, therefore won't affect the batch corrections but
#' adjust the results in down-stream analyses.
#' @param sv.number number of surrogate variables to be calculated using SVA, best left untouched.
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
#' @param verbose prints messages through running the pipeline
#'
#' @examples
#' \donttest{
#' data(atac_seq_consensus_bm) # calls 'bed'
#'
#' # a vector for comparing the examples
#' contrasts <- sapply(strsplit(colnames(bed), split = "-", fixed = TRUE),
#'                     function(x){x[1]})[4:25]
#'
#' results <- cinaR(bed, contrasts, reference.genome = "mm10")
#' }
#'
#'
#' @return returns differentially accessible peaks
#'
#' @export
cinaR <-
  function(matrix,
           contrasts,
           experiment.type = "ATAC-Seq",
           DA.choice = 1,
           DA.fdr.threshold = 0.05,
           DA.lfc.threshold = 0,
           comparison.scheme = "OVO",
           save.DA.peaks = FALSE,
           DA.peaks.path = NULL,
           norm.method = "cpm",
           filter.method = "custom",
           library.threshold = 2,
           cpm.threshold = 1,
           TSS.threshold = 50e3,
           show.annotation.pie = FALSE,
           reference.genome = NULL,
           batch.correction = FALSE,
           batch.information = NULL,
           additional.covariates = NULL,
           sv.number = NULL,
           run.enrichment = TRUE,
           enrichment.method = NULL,
           enrichment.FDR.cutoff = 1,
           background.genes.size = 20e3,
           geneset = NULL,
           verbose = TRUE) {

    # Printing function
    verbosePrint <- verboseFn(verbose)

    verbosePrint(">> Experiment type: ", experiment.type)


    if (is.null(reference.genome)) {
      warning("'reference.genome' is not set, therefore hg38 will be used!")
      reference.genome <- "hg38"
    }

    if (experiment.type == "ATAC-Seq"){

      if (length(contrasts) != (ncol(matrix) - 3)) {
        stop("Length of 'contrasts' must be equal to number of samples in 'matrix'")
      }

      # collapse chr, start, end and make them rownames
      cp.rownames <-
        apply(matrix[, 1:3], 1, function(x) {
          paste0(trimws(x), collapse = "_")
        })
      cp <- matrix[, -c(1:3)]
      rownames(cp) <- cp.rownames

      # filter low expressed peaks
      cp.filtered <-
        filterConsensus(cp, library.threshold = library.threshold, cpm.threshold = cpm.threshold)

      verbosePrint(">> Matrix is filtered!")

      # annotate the peaks to the closest TSS
      cp.filtered.annotated <- annotatePeaks(cp.filtered,
                                             reference.genome = reference.genome,
                                             show.annotation.pie = show.annotation.pie,
                                             verbose = verbose)

      # filter distance to TSS
      final.matrix <-
        cp.filtered.annotated[abs(cp.filtered.annotated$distanceToTSS) <= TSS.threshold, ]


    } else if (experiment.type == "RNA-Seq"){

      if (length(contrasts) != (ncol(matrix) - 1)) {
        stop("Length of 'contrasts' must be equal to number of samples in `matrix`")
      }

      verbosePrint(">> Arranging count matrix...")

      # Remove spike-ins for now (may not exists in your data)
      matrix <- matrix [!grepl("^ERCC", matrix[,1]),]

      # Eliminate any homologs
      # TODO could be dangerous to do it this way, find a better version...
      matrix[,1] <- sapply(strsplit(matrix[,1], ".", fixed = TRUE), function(x){x[1]})

      # Order genes according to their standard deviation in decreasing order
      matrix <- matrix [rev(order(apply(matrix[,-1], 1, stats::sd))),]

      # Remove duplicated genes
      matrix <- matrix [!duplicated(matrix[,1]),]

      # Make the gene names the row names
      rownames(matrix) <- matrix[,1]

      # Filter the genes
      matrix <- matrix[,-1]

      # Enforce all counts to be integers
      matrix <- round(matrix, 0)

      # filter low expressed peaks
      final.matrix <-
        filterConsensus(matrix, library.threshold = library.threshold, cpm.threshold = cpm.threshold)

      verbosePrint(">> Matrix is filtered!")

    } else {
      stop("`experiment.type` must be either 'ATAC-Seq' or 'RNA-Seq'")
    }


    if (!is.null(enrichment.method)) {
      if (enrichment.method == "GSEA" & run.enrichment == TRUE) {
        verbosePrint(
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
        final.matrix = final.matrix,
        contrasts = contrasts,
        experiment.type = experiment.type,
        DA.choice = DA.choice,
        DA.fdr.threshold = DA.fdr.threshold,
        DA.lfc.threshold = DA.lfc.threshold,
        comparison.scheme = comparison.scheme,
        save.DA.peaks = save.DA.peaks,
        DA.peaks.path = DA.peaks.path,
        batch.correction = batch.correction,
        batch.information = batch.information,
        additional.covariates = additional.covariates,
        sv.number = sv.number,
        verbose = verbose
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
          experiment.type = experiment.type,
          reference.genome = reference.genome,
          enrichment.method = enrichment.method,
          enrichment.FDR.cutoff = enrichment.FDR.cutoff,
          background.genes.size = background.genes.size,
          verbose = verbose
        )
      verbosePrint(">> Enrichment results are ready...")
      verbosePrint(">> Done!")
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
#' @param cp consensus peak matrix, with unique ids at rownames.
#' @param filter.method filtering method for low expressed peaks
#' @param library.threshold number of libraries a peak occurs so that it is not filtered default set to 2
#' @param cpm.threshold count per million threshold for not to filter a peak
#'
#' @return returns differentially accessible peaks
#'
#' @examples
#' set.seed(123)
#' cp <- matrix(rexp(200, rate=.1), ncol=20)
#'
#' ## using cpm function from `edgeR` package
#' cp.filtered <- filterConsensus(cp)
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
    return(cp.filtered)
  }


#' normalizeConsensus
#'
#' Normalizes consensus peak using different methods
#'
#' @param cp bed formatted consensus peak matrix: CHR, START, STOP and raw peak counts (peaks by 3+samples)
#' @param norm.method normalization method for consensus peaks
#' @param log.option logical, log option for cpm function in edgeR
#' @return Normalized consensus peaks
#' @examples
#'
#' set.seed(123)
#' cp <- matrix(rexp(200, rate=.1), ncol=20)
#'
#' ## using cpm function from `edgeR` package
#' cp.normalized <- normalizeConsensus(cp)
#'
#' ## quantile normalization option
#' cp.normalized <- normalizeConsensus(cp, norm.method = "quantile")
#' @export
normalizeConsensus <-
  function(cp, norm.method = "cpm", log.option = FALSE) {
    if (norm.method == "cpm") {
      if (log.option){
        # we don't use cpm log option,
        # to make sure it does not yield any negative values.
        cp.norm <- log2(edgeR::cpm(cp) + 1)
      } else {
        cp.norm <- edgeR::cpm(cp)
      }

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
#' @param verbose prints messages through running the pipeline
#'
#' @return DApeaks returns DA peaks
annotatePeaks <-
  function(cp,
           reference.genome,
           show.annotation.pie = FALSE,
           verbose) {
    

    bed <-
      as.data.frame(do.call(rbind, strsplit(rownames(cp), "_", fixed = TRUE)))
    colnames(bed) <- c("CHR", "Start", "End")
    bed.GRanges <- GenomicRanges::GRanges(bed)

    if (reference.genome == "hg38") {
      if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
        
        message(
          "Package \"TxDb.Hsapiens.UCSC.hg38.knownGene\" needed for this
             function to work. Please install it."
        )
        return(NULL)
        
      }
      txdb <-
        TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
      genome <- cinaR::grch38
      reference.genome <- "hg38"
    } else if (reference.genome == "hg19") {
      if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE)) {
        message(
          "Package \"TxDb.Hsapiens.UCSC.hg19.knownGene\" needed for this
             function to work. Please install it."
        )
        return(NULL)
      }
      txdb <-
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
      genome <- cinaR::grch37
    } else if (reference.genome == "mm10") {
      if (!requireNamespace("TxDb.Mmusculus.UCSC.mm10.knownGene", quietly = TRUE)) {
        message(
          "Package \"TxDb.Mmusculus.UCSC.mm10.knownGene\" needed for this
             function to work. Please install it.",
          call. = FALSE
        )
        return(NULL)
      }
      txdb <-
        TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
      genome <- cinaR::grcm38
    } else {
      stop ("reference.genome should be 'hg38', 'hg19' or 'mm10'")
    }

    # annotate peaks
    annoPeaks <- ChIPseeker::annotatePeak(bed.GRanges, TxDb = txdb, verbose = verbose)


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

    return(cbind(data.frame(annoPeaks.anno), cp))
  }

#' Differential Analyses
#'
#' Runs differential analyses pipeline of choice on consensus peaks
#'
#' @param final.matrix Annotated Consensus peaks
#' @param contrasts user-defined contrasts for comparing samples
#' @param experiment.type The type of experiment either set to "ATAC-Seq" or "RNA-Seq"
#' @param DA.choice determines which pipeline to run:
#' (1) edgeR, (2) limma-voom, (3) limma-trend, (4) DEseq2
#' @param DA.fdr.threshold fdr cut-off for differential analyses
#' @param DA.lfc.threshold log-fold change cutoff for differential analyses
#' @param comparison.scheme either one-vs-one (OVO) or one-vs-all (OVA) comparisons.
#' @param save.DA.peaks logical, saves differentially accessible peaks to an excel file
#' @param DA.peaks.path the path which the excel file of the DA peaks will be saved,
#' if not set it will be saved to current directory.
#' @param batch.correction logical, if set will run unsupervised batch correction
#' via sva (default) or if the batch information is known `batch.information`
#' argument should be provided by user.
#' @param batch.information character vector, given by user.
#' @param additional.covariates vector or data.frame, this parameter will be directly added to design
#' matrix before running the differential analyses, therefore won't affect the batch corrections but
#' adjust the results in down-stream analyses.
#' @param sv.number number of surrogate variables to be calculated using SVA, best left untouched.
#' @param verbose prints messages through running the pipeline
#' @return returns consensus peaks (batch corrected version if enabled) and DA peaks
differentialAnalyses <- function(final.matrix,
                                 contrasts,
                                 experiment.type,
                                 DA.choice,
                                 DA.fdr.threshold,
                                 DA.lfc.threshold,
                                 comparison.scheme,
                                 save.DA.peaks,
                                 DA.peaks.path,
                                 batch.correction,
                                 batch.information,
                                 additional.covariates,
                                 sv.number,
                                 verbose) {

  # Printing function
  verbosePrint <- verboseFn(verbose)

  # silence CRAN build notes
  log2FoldChange <- padj <- NULL

  if (experiment.type == "ATAC-Seq"){

    cp.meta <- final.matrix[, 1:15]
    cp.metaless <- final.matrix[, 16:ncol(final.matrix)]

  } else { # RNA-seq
    cp.metaless <- final.matrix
  }

  design <- stats::model.matrix(~ 0 + contrasts)

  # https://www.biostars.org/p/461026/
  if (batch.correction) {

    if (is.null(batch.information)) {
      ## First normalize the consensus peaks to avoid detecting the effects
      ## confounding from library size as Michael Love and Jeff Leek suggests
      ## in this thread:
      verbosePrint(">> Running SVA for batch correction...")

      cp.metaless.normalized <- normalizeConsensus(cp.metaless, log.option = TRUE)
      mod  <- stats::model.matrix(~ 0 + contrasts)
      mod0 <- cbind(rep(1, length(contrasts)))

      # calculate the batch effects
      if (is.null(sv.number)){
        sva.res <-
          sva::svaseq(cp.metaless.normalized, mod, mod0)
      } else {
        sva.res <-
          sva::svaseq(cp.metaless.normalized, mod, mod0, n.sv = sv.number)
      }

      # batch effect additional matrix
      add.batch <- sva.res$sv

      # make the colnames prettier just for fun
      colnames(add.batch) <- paste0("SV", c(1:ncol(add.batch)))

      # add it to the design matrix
      design <-
        cbind(design, add.batch)

      # batch corrected/normalized consensus peaks created for PCA/Heatmaps
      cp.batch.corrected <- limma::removeBatchEffect(cp.metaless.normalized, covariates = add.batch)

    } else {
      # if there is batch information available
      verbosePrint(">> Adding batch information to design matrix...")

      if (nrow(design) != length(batch.information)) {
        stop("Number of samples and `batch.information` should be same length!")
      }

      cp.metaless.normalized <- normalizeConsensus(cp.metaless, log.option = TRUE)

      design <- cbind(design, BatchInfo = batch.information)

      # batch corrected consensus peaks created for PCA/Heatmaps
      cp.batch.corrected <- limma::removeBatchEffect(cp.metaless.normalized, batch = batch.information)
    }
  }

  if (!is.null(additional.covariates)) {

    # If additional covariates are already data.frame
    # this line does not change anything!
    df.covariates <- data.frame(additional.covariates)

    if (nrow(df.covariates) != nrow(design)){
      stop("Number of samples in `additional.covariates` should match the with the sample size!")
    }
    design <- cbind(design, additional.covariates)

    verbosePrint(">> Additional covariates are added to design matrix...")
  }

  # Add intercept term for multiple comparisons

  rownames(design) <- colnames(cp.metaless)
  colnames(design) <- gsub("contrasts", "", colnames(design))


  if (comparison.scheme == "OVO"){ # one vs one

    # Create contrasts for all comparisons
    combs <-
      utils::combn(colnames(design)[1:length(unique(contrasts))], 2)

    contrasts.order <- c(1:length(unique(contrasts)))
    names(contrasts.order) <- unique(contrasts)

    # Re-order contrasts according to group order
    combs <- apply(combs, 2, function(x){
      names(sort(contrasts.order[x]))
    })

    cc <- apply(combs, 2,
                function(x) {
                  paste0(paste(x, collapse = "_"), "=", x[1], "-", x[2])
                })

  } else if(comparison.scheme == "OVA") { # one vs all
    comps <- colnames(design)[1:length(unique(contrasts))]
    cc <- NULL
    for (i in seq(1,length(comps))){
      cc <- cbind(cc,
                  paste0(comps[i],"_REST=",
                         comps[i], "-",
                         paste0(comps[-i], "/", length(comps)-1, collapse = "-")))
    }
  } else {
    stop("Comparison scheme should be either 'OVO' (one vs one) or 'OVA' (one vs all)")
  }

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

    verbosePrint(
      ">> Method: edgeR\n\tFDR:",
      DA.fdr.threshold,
      "& abs(logFC)<",
      DA.lfc.threshold
    )

    y <- edgeR::DGEList(counts = cp.metaless, group = contrasts)

    # Calculate normalization factors for library sizes with TMM
    y <- edgeR::calcNormFactors(y, method = "TMM")

    # Estimate dispersion for genes with Bayesian Shrinkage
    verbosePrint(">> Estimating dispersion...")
    y <- edgeR::estimateDisp(y, design)

    # Fit the model
    verbosePrint(">> Fitting GLM...")
    fit.glm <- edgeR::glmQLFit(y, design)

    for (i in seq_len(ncol(ccc))) {
      contrast.name <- colnames(ccc)[i]
      qlf <- edgeR::glmQLFTest(fit.glm, contrast = ccc[, i])
      # plotMD(qlf, main = contrast.name, p.value = 0.1)
      top.table <-
        edgeR::topTags(qlf, n = Inf, p.value = DA.fdr.threshold)$table

      # ifelse does not return the dataframe for some reason,
      # therefore, implemented this check explicitly
      if(is.null(top.table)){
        top.table <- data.frame()
      }

      if (nrow(top.table) > 0) {

        top.table <- top.table[abs(top.table$logFC) >= DA.lfc.threshold,]

        if(experiment.type == "ATAC-Seq"){
          top.table <- merge(cp.meta, top.table, by = 0)
          # Refactor to uniformize DA results
          top.table <- top.table[, c(1:17, 21)]
        } else {
          top.table <- top.table[, c(1, 5)]
          top.table <- cbind(gene_name = rownames(top.table), top.table)
          rownames(top.table) <- NULL
        }

        DA.peaks[[contrast.name]] <- top.table

      } else {
        DA.peaks[[contrast.name]] <- list()
      }
    }
  } else if (DA.choice == 2) {
    ## limma-voom
    verbosePrint(
      ">> Method: limma-voom\n\tFDR:",
      DA.fdr.threshold,
      "& abs(logFC)<",
      DA.lfc.threshold
    )
    v <- limma::voom(cp.metaless, design, plot = FALSE)
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

      if(experiment.type == "ATAC-Seq"){
        top.table <- merge(cp.meta, top.table, by = 0)

        # Refactor to uniform DA results
        top.table <- top.table[, c(1:17, 21)]
        colnames(top.table)[18] <- "FDR"
      } else {
        # Refactor to uniform DA results
        top.table <- top.table[, c(1, 5)]
        colnames(top.table)[2] <- "FDR"
        top.table<- cbind(gene_name = rownames(top.table), top.table)
      }

      # Safety check
      if(is.null(top.table)){
        top.table <- data.frame()
      }

      if (nrow(top.table) > 0) {
        rownames(top.table) <- NULL
        DA.peaks[[contrast.name]] <- top.table
      } else {
        DA.peaks[[contrast.name]] <- list()
      }

    }
  } else if (DA.choice == 3) {
    ## limma-trend
    verbosePrint(
      ">> Method: limma-trend\n\tFDR:",
      DA.fdr.threshold,
      "& abs(logFC)<",
      DA.lfc.threshold
    )
    fit.trend <- limma::lmFit(cp.metaless, design)
    fit.trend2 <-
      limma::eBayes(limma::contrasts.fit(fit.trend, ccc),
                    trend = TRUE)
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

      if (experiment.type == "ATAC-Seq"){
        top.table <- merge(cp.meta, top.table, by = 0)

        # Refactor to uniformize DA results
        top.table <- top.table[, c(1:17, 21)]
        colnames(top.table)[18] <- "FDR"
      } else {
        # Refactor to uniform DA results
        top.table <- top.table[, c(1, 5)]
        colnames(top.table)[2] <- "FDR"
        top.table<- cbind(gene_name = rownames(top.table), top.table)
      }

      # Safety check
      if(is.null(top.table)){
        top.table <- data.frame()
      }
      if (nrow(top.table) > 0) {
        rownames(top.table) <- NULL
        DA.peaks[[contrast.name]] <- top.table
      } else {
        DA.peaks[[contrast.name]] <- list()
      }
    }
  } else if (DA.choice == 4) {
    ## DEseq2
    verbosePrint(
      ">> Method: DEseq2\n\tFDR:",
      DA.fdr.threshold,
      "& abs(logFC)<",
      DA.lfc.threshold
    )


    # Assign each sample to its group
    colData <- as.data.frame(cbind(colnames(cp.metaless), contrasts))
    colnames(colData)  = c("sample", "groups")

    # Create DEseq Object
    dds <-
      DESeq2::DESeqDataSetFromMatrix(countData = cp.metaless,
                                     colData = colData,
                                     design = ~ groups)

    dds = DESeq2::DESeq(dds, parallel = TRUE)

    # Create DE gene list for DESeq2

    for (i in seq_len(ncol(ccc))) {
      contrast.name <- colnames(ccc)[i]
      DEseq.contrast <- rownames(ccc)[ccc[, i] != 0]
      res <-
        DESeq2::results(
          dds,
          c("groups", DEseq.contrast[2], DEseq.contrast[1]),
          parallel = TRUE,
          tidy = TRUE
        )
      rownames(res) <- res$row
      res.ordered <- res[order(res$pvalue), ]
      res.significant <-
        subset(res.ordered,
               padj <= DA.fdr.threshold &
                 abs(log2FoldChange) >= DA.lfc.threshold)

      if (experiment.type == "ATAC-Seq"){
        res.significant <- merge(cp.meta, res.significant, by = 0)
        top.table <- res.significant[, c(1:16, 19, 23)]
        colnames(top.table)[c(17, 18)] <- c("logFC", "FDR")
      } else {
        top.table <- res.significant[, c(1,3,7)]
        colnames(top.table) <- c("gene_name", "logFC", "FDR")
      }

      if(is.null(top.table)){
        top.table <- data.frame()
      }

      if (nrow(top.table) > 0) {
        rownames(top.table) <- NULL
        DA.peaks[[contrast.name]] <- top.table
      } else {
        DA.peaks[[contrast.name]] <- list()
      }
    }
  }

  verbosePrint(">> DA peaks are found!")

  if (save.DA.peaks) {
    # Make sure every list
    DA.peaks.dfs <- lapply(DA.peaks, data.frame)

    if (is.null(DA.peaks.path)) {
      verbosePrint(">> Saving DA peaks to current directory as DApeaks.xlsx...")
      writexl::write_xlsx(x = DA.peaks.dfs, path = "./DApeaks.xlsx")
    } else {
      verbosePrint(paste0(">> Saving DA peaks to ", DA.peaks.path, "..."))
      writexl::write_xlsx(x = DA.peaks.dfs, path = DA.peaks.path)
    }
  }

  if (batch.correction){
    return(list (cp = cp.batch.corrected, DA.peaks = DA.peaks))
  }
  return(list (cp = normalizeConsensus(cp.metaless, log.option = T), DA.peaks = DA.peaks))
}
