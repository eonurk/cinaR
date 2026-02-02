# normalizeConsensus

Normalizes consensus peak using different methods

## Usage

``` r
normalizeConsensus(cp, norm.method = "cpm", log.option = FALSE)
```

## Arguments

- cp:

  bed formatted consensus peak matrix: CHR, START, STOP and raw peak
  counts (peaks by 3+samples)

- norm.method:

  normalization method for consensus peaks

- log.option:

  logical, log option for cpm function in edgeR

## Value

Normalized consensus peaks

## Examples

``` r
set.seed(123)
cp <- matrix(rexp(200, rate=.1), ncol=20)

## using cpm function from `edgeR` package
cp.normalized <- normalizeConsensus(cp)

## quantile normalization option
cp.normalized <- normalizeConsensus(cp, norm.method = "quantile")
```
