if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("maftools", force=TRUE)

library(maftools)

query <- GDCquery(
  project = "TCGA-STAD", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  legacy = FALSE, 
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query)
maf <- GDCprepare(query)

maf <- maf %>% maftools::read.maf

maf<-read.maf(maf)

datatable(getSampleSummary(maf),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)


oncoplot(maf = maf, top = 10, removeNonMutated = TRUE)
titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = titv)


