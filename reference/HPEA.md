# HPEA

Hyper-geometric p-value enrichment analyses, looking for
over-representation of a set of genes on given pathways.

## Usage

``` r
HPEA(genes, geneset, background.genes.size)
```

## Arguments

- genes:

  DA gene names to be checked if they are over-represented or not.

- geneset:

  Pathways to be used in enrichment analyses. If not set vp2008
  (Chaussabel, 2008) immune modules will be used. This can be set to any
  geneset using \`read.gmt\` function from \`qusage\` package. Different
  modules are available: https://www.gsea-msigdb.org/gsea/downloads.jsp.

- background.genes.size:

  number of background genes for hyper-geometric p-value calculations.
  Default is 20,000.

## Value

data.frame, list of pathways and their enrichment (adjusted) p-values.

## Examples

``` r
# \donttest{
library(cinaR)

data("VP2008")
genes.to.test <- vp2008[[1]][1:10]
HPEA(genes.to.test,vp2008, background.genes.size = 20e3)
#>                                                       module.name        p.val
#> Plasma cells                                         Plasma cells 1.853276e-26
#> Erythrocytes                                         Erythrocytes 1.251584e-04
#> Platelets                                               Platelets 1.000000e+00
#> B cells                                                   B cells 1.000000e+00
#> U_cAMP/NF-KB activation                   U_cAMP/NF-KB activation 1.000000e+00
#> Myeloid lineage 1                               Myeloid lineage 1 1.000000e+00
#> U_P53 signaling                                   U_P53 signaling 1.000000e+00
#> MHC/Ribosomal proteins                     MHC/Ribosomal proteins 1.000000e+00
#> U_metabolism/replication                 U_metabolism/replication 1.000000e+00
#> Cytotoxic cells                                   Cytotoxic cells 1.000000e+00
#> Neutrophils                                           Neutrophils 1.000000e+00
#> Ribosomal proteins                             Ribosomal proteins 1.000000e+00
#> U_Immunity/cytoskeleton                   U_Immunity/cytoskeleton 1.000000e+00
#> Myeloid lineage 2                               Myeloid lineage 2 1.000000e+00
#> Unknown                                                   Unknown 1.000000e+00
#> T Cells                                                   T Cells 1.000000e+00
#> U_T cells/cytoskeleton                     U_T cells/cytoskeleton 1.000000e+00
#> U_Immsurface/cytokines/signaling U_Immsurface/cytokines/signaling 1.000000e+00
#> U_RAS/kinases                                       U_RAS/kinases 1.000000e+00
#> Interferon-inducible                         Interferon-inducible 1.000000e+00
#> Inflammation I                                     Inflammation I 1.000000e+00
#> Inflammation II                                   Inflammation II 1.000000e+00
#> U_protphosphatases/PI3K                   U_protphosphatases/PI3K 1.000000e+00
#> U_hemoglobin                                         U_hemoglobin 1.000000e+00
#> U_mitochondrial proteins                 U_mitochondrial proteins 1.000000e+00
#> U_proteasome/ubiquitin cx               U_proteasome/ubiquitin cx 1.000000e+00
#> U_enzymes                                               U_enzymes 1.000000e+00
#> U_kinases/phosphatases                     U_kinases/phosphatases 1.000000e+00
#>                                                                                   overlapping.genes
#> Plasma cells                     IGHD,PLEKHG1,C13orf18,ZDHHC23,HLA-DOA,MS4A1,TCF4,HTPAP,SMC6L1,IGHM
#> Erythrocytes                                                                              IGHD,IGHM
#> Platelets                                                                                          
#> B cells                                                                                            
#> U_cAMP/NF-KB activation                                                                            
#> Myeloid lineage 1                                                                                  
#> U_P53 signaling                                                                                    
#> MHC/Ribosomal proteins                                                                             
#> U_metabolism/replication                                                                           
#> Cytotoxic cells                                                                                    
#> Neutrophils                                                                                        
#> Ribosomal proteins                                                                                 
#> U_Immunity/cytoskeleton                                                                            
#> Myeloid lineage 2                                                                                  
#> Unknown                                                                                            
#> T Cells                                                                                            
#> U_T cells/cytoskeleton                                                                             
#> U_Immsurface/cytokines/signaling                                                                   
#> U_RAS/kinases                                                                                      
#> Interferon-inducible                                                                               
#> Inflammation I                                                                                     
#> Inflammation II                                                                                    
#> U_protphosphatases/PI3K                                                                            
#> U_hemoglobin                                                                                       
#> U_mitochondrial proteins                                                                           
#> U_proteasome/ubiquitin cx                                                                          
#> U_enzymes                                                                                          
#> U_kinases/phosphatases                                                                             
#>                                  module.overlapping.ratio        adj.p
#> Plasma cells                                   0.17241379 5.189173e-25
#> Erythrocytes                                   0.05882353 1.752218e-03
#> Platelets                                      0.00000000 1.000000e+00
#> B cells                                        0.00000000 1.000000e+00
#> U_cAMP/NF-KB activation                        0.00000000 1.000000e+00
#> Myeloid lineage 1                              0.00000000 1.000000e+00
#> U_P53 signaling                                0.00000000 1.000000e+00
#> MHC/Ribosomal proteins                         0.00000000 1.000000e+00
#> U_metabolism/replication                       0.00000000 1.000000e+00
#> Cytotoxic cells                                0.00000000 1.000000e+00
#> Neutrophils                                    0.00000000 1.000000e+00
#> Ribosomal proteins                             0.00000000 1.000000e+00
#> U_Immunity/cytoskeleton                        0.00000000 1.000000e+00
#> Myeloid lineage 2                              0.00000000 1.000000e+00
#> Unknown                                        0.00000000 1.000000e+00
#> T Cells                                        0.00000000 1.000000e+00
#> U_T cells/cytoskeleton                         0.00000000 1.000000e+00
#> U_Immsurface/cytokines/signaling               0.00000000 1.000000e+00
#> U_RAS/kinases                                  0.00000000 1.000000e+00
#> Interferon-inducible                           0.00000000 1.000000e+00
#> Inflammation I                                 0.00000000 1.000000e+00
#> Inflammation II                                0.00000000 1.000000e+00
#> U_protphosphatases/PI3K                        0.00000000 1.000000e+00
#> U_hemoglobin                                   0.00000000 1.000000e+00
#> U_mitochondrial proteins                       0.00000000 1.000000e+00
#> U_proteasome/ubiquitin cx                      0.00000000 1.000000e+00
#> U_enzymes                                      0.00000000 1.000000e+00
#> U_kinases/phosphatases                         0.00000000 1.000000e+00
# }
```
