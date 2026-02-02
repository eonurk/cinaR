# annotatePeaks

Runs DA pipeline and makes it ready for enrichment analyses

## Usage

``` r
annotatePeaks(cp, reference.genome, show.annotation.pie = FALSE, verbose)
```

## Arguments

- cp:

  bed formatted consensus peak matrix: CHR, START, STOP and raw peak
  counts (peaks by 3+samples)

- reference.genome:

  genome of interested species. It should be 'hg38', 'hg19' or 'mm10'.

- show.annotation.pie:

  shows the annotation pie chart produced with ChipSeeker

- verbose:

  prints messages through running the pipeline

## Value

DApeaks returns DA peaks
