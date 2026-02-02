# filterConsensus

Filters lowly expressed peaks from down-stream analyses

## Usage

``` r
filterConsensus(
  cp,
  filter.method = "custom",
  library.threshold = 2,
  cpm.threshold = 1
)
```

## Arguments

- cp:

  consensus peak matrix, with unique ids at rownames.

- filter.method:

  filtering method for low expressed peaks

- library.threshold:

  number of libraries a peak occurs so that it is not filtered default
  set to 2

- cpm.threshold:

  count per million threshold for not to filter a peak

## Value

returns differentially accessible peaks

## Examples

``` r
set.seed(123)
cp <- matrix(rexp(200, rate=.1), ncol=20)

## using cpm function from `edgeR` package
cp.filtered <- filterConsensus(cp)
```
