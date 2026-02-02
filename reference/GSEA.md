# GSEA

Gene set enrichment analyses, runs 'fgsea' package implementation with
preset values.

## Usage

``` r
GSEA(genes, geneset)
```

## Arguments

- genes:

  DA gene names to be checked if they are over-represented or not.

- geneset:

  Pathways to be used in enrichment analyses. If not set vp2008
  (Chaussabel, 2008) immune modules will be used. This can be set to any
  geneset using \`read.gmt\` function from \`qusage\` package. Different
  modules are available: https://www.gsea-msigdb.org/gsea/downloads.jsp.

## Value

data.frame, list of pathways and their enrichment (adjusted) p-values.

## References

G. Korotkevich, V. Sukhov, A. Sergushichev. Fast gene set enrichment
analysis. bioRxiv (2019), doi:10.1101/060012

## Examples

``` r
# \donttest{
library(cinaR)
library(fgsea)
data(examplePathways)
data(exampleRanks)
GSEA(exampleRanks, examplePathways)
#>                                                                                      pathway
#>                                                                                       <char>
#>   1:                                                                1221633_Meiotic_Synapsis
#>   2:                                   1445146_Translocation_of_Glut4_to_the_Plasma_Membrane
#>   3: 442533_Transcriptional_Regulation_of_Adipocyte_Differentiation_in_3T3-L1_Pre-adipocytes
#>   4:                                                                  508751_Circadian_Clock
#>   5:                                               5334727_Mus_musculus_biological_processes
#>  ---                                                                                        
#> 582:                                                          6096057_Adaptive_Immune_System
#> 583:                          6096937_Class_I_MHC_mediated_antigen_processing_&_presentation
#> 584:                                                            6096957_ER-Phagosome_pathway
#> 585:                                           6096958_Antigen_processing-Cross_presentation
#> 586:                                                            912497_Meiotic_Recombination
#>           pval       padj    log2err         ES        NES  size
#>          <num>      <num>      <num>      <num>      <num> <int>
#>   1: 0.5539305 0.73273881 0.06928365  0.2885754  0.9418916    27
#>   2: 0.6684211 0.81944506 0.05822162  0.2387284  0.8590286    39
#>   3: 0.1192661 0.28410532 0.19578900 -0.3640706 -1.3024359    31
#>   4: 0.8044280 0.90478855 0.05195125  0.2516324  0.7391437    17
#>   5: 0.3813953 0.58694683 0.07998588  0.2469065  1.0453809   106
#>  ---                                                            
#> 582: 0.7088608 0.84142305 0.05688642  0.2568888  0.8256794    25
#> 583: 0.7088608 0.84142305 0.05688642  0.2568888  0.8256794    25
#> 584: 0.3613139 0.56986701 0.09196861  0.3477495  1.1006039    23
#> 585: 0.7088608 0.84142305 0.05688642  0.2568888  0.8256794    25
#> 586: 0.0248887 0.09349217 0.35248786  0.5722386  1.7116914    18
#>                                                                               leadingEdge
#>                                                                                    <char>
#>   1:                                                                    15270,12189,71846
#>   2:                          17918,19341,20336,22628,22627,20619,16579,16568,11651,12315
#>   3:            76199,19014,26896,229003,17977,17978,12537,70208,67381,59024,327987,20602
#>   4:                                                                          20893,59027
#>   5:             60406,19361,15270,20893,12189,68240,71846,20018,192191,12567,19891,59027
#>  ---                                                                                     
#> 582: 26440,26444,26445,19170,26442,19177,667803,12317,53421,19172,19166,19173,16913,26441
#> 583: 26440,26444,26445,19170,26442,19177,667803,12317,53421,19172,19166,19173,16913,26441
#> 584: 26440,26444,26445,19170,26442,19177,667803,12317,53421,19172,19166,19173,16913,26441
#> 585: 26440,26444,26445,19170,26442,19177,667803,12317,53421,19172,19166,19173,16913,26441
#> 586:                                                  19361,15270,12189,68240,12567,19891
# }
```
