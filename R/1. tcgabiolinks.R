if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
install.packages("DT")

library(TCGAbiolinks)
library(dplyr)
library(DT)

query <- GDCquery(
  project = c("TCGA-GBM", "TCGA-LGG"),
  data.category = "DNA Methylation",
  legacy = FALSE,
  platform = c("Illumina Human Methylation 450"),
  sample.type = "Recurrent Tumor"
)
datatable(
  getResults(query), 
  filter = 'top',
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = FALSE
)


subtypes <- PanCancerAtlas_subtypes()
DT::datatable(subtypes,
              filter = 'top',
              options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
              rownames = FALSE)



COAD.subtype <- TCGAquery_subtype(tumor = c("coad"))
stad.subtype <- TCGAquery_subtype(tumor = c("stad"))
read.subtype <- TCGAquery_subtype(tumor = c("read"))


View(subtypes)

query <- GDCquery(
  project = "TCGA-CHOL", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  legacy = FALSE, 
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query)
maf <- GDCprepare(query)

datatable(maf[1:20,],
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

