# Differential Analyses

Runs differential analyses pipeline of choice on consensus peaks

## Usage

``` r
differentialAnalyses(
  final.matrix,
  contrasts,
  experiment.type,
  DA.choice,
  DA.fdr.threshold,
  DA.lfc.threshold,
  comparison.scheme,
  save.DA.peaks,
  DA.peaks.path,
  norm.method,
  batch.correction,
  batch.information,
  additional.covariates,
  sv.number,
  verbose
)
```

## Arguments

- final.matrix:

  Annotated Consensus peaks

- contrasts:

  user-defined contrasts for comparing samples

- experiment.type:

  The type of experiment either set to "ATAC-Seq" or "RNA-Seq"

- DA.choice:

  determines which pipeline to run: (1) edgeR, (2) limma-voom, (3)
  limma-trend, (4) DEseq2

- DA.fdr.threshold:

  fdr cut-off for differential analyses

- DA.lfc.threshold:

  log-fold change cutoff for differential analyses

- comparison.scheme:

  either one-vs-one (OVO) or one-vs-all (OVA) comparisons.

- save.DA.peaks:

  logical, saves differentially accessible peaks to an excel file

- DA.peaks.path:

  the path which the excel file of the DA peaks will be saved, if not
  set it will be saved to current directory.

- norm.method:

  normalization method for consensus peaks

- batch.correction:

  logical, if set will run unsupervised batch correction via sva
  (default) or if the batch information is known \`batch.information\`
  argument should be provided by user.

- batch.information:

  character vector, given by user.

- additional.covariates:

  vector or data.frame, this parameter will be directly added to design
  matrix before running the differential analyses, therefore won't
  affect the batch corrections but adjust the results in down-stream
  analyses.

- sv.number:

  number of surrogate variables to be calculated using SVA, best left
  untouched.

- verbose:

  prints messages through running the pipeline

## Value

returns consensus peaks (batch corrected version if enabled) and DA
peaks
